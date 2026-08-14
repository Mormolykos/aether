# Aether

A real-time aerospace telemetry engine in Rust: it ingests live ADS-B aircraft
broadcasts, maintains a filtered track picture, and screens every pair of tracks for
closest approach against ICAO separation minima.

Everything below is separated into **measured** results — numbers produced by running
this code against the live feed — and **architectural** descriptions of how it works.
Where something has not been measured, it says so.

---

## What Aether does

- Polls a public ADS-B feed asynchronously and turns each response into strongly typed
  contacts, discarding malformed and incomplete records at the boundary.
- Places each measurement at the moment it was **observed**, not the moment it arrived.
- Converts WGS84 geodetic positions into a local East/North/Up frame in metres.
- Maintains one track per ICAO address, estimating position and velocity with a
  constant-velocity Kalman filter per axis.
- Rejects measurements that are kinematically implausible, and counts the rejections.
- Ages tracks out on the staleness of their evidence.
- Screens all firm pairs for closest point of approach, in closed form.
- Exposes the closest-approach primitive over a C ABI for use from C or C++.
- Bounds the memory a single upstream response can cause it to allocate.

It does **not** contain targeting, engagement, weapon, or fire-control functionality of
any kind, and is not intended to.

## Running it

```sh
cp .env.example .env
cargo run --release
```

Configuration is read from `.env` in the working directory, or from a path given as the
first argument. It holds the observer position, poll and cycle periods, filter tuning,
separation minima, and the response-body limit.

**Aether needs no credentials.** The ADS-B feed requires no API key, and there is nothing
secret in `.env.example`. The real `.env` is still git-ignored — a local file is where a
different observer position, an alternative feed, or a tuning experiment ends up, and
none of those belong in a commit.

Press Ctrl-C to stop; it prints a session summary on the way out.

---

## Architecture

```text
 sensor task                          tracker task
┌──────────────────┐   bounded    ┌────────────────────────────────────┐
│ poll feed        │   mpsc(8)    │ associate by ICAO                  │
│ bound body size  │ ──Batch───▶  │ place measurement at observed time │
│ decode → Contact │              │ predict → gate → update            │
│ count health     │              │ age out stale tracks               │
└──────────────────┘              │ screen pairs → render              │
                                  └────────────────────────────────────┘
```

Two Tokio tasks, one channel, no shared mutable state. The sensor task owns its HTTP
client and its counters; the tracker task owns every track. Nothing is behind a mutex,
so the fixed-rate loop has no lock to contend on and a slow network read cannot stall
the picture.

Sensor health is propagated through the existing channel by riding on each batch, rather
than through shared state — there is no `Arc<Mutex<_>>` on this path. That is a
statement about this design, not a claim that the system is lock-free in the
formal, progress-guarantee sense.

The channel is bounded deliberately. When the tracker falls behind, the sensor sheds a
batch instead of growing an unbounded queue of increasingly stale observations.

## Telemetry ingestion

ADS-B is a cooperative surveillance broadcast: transponder-equipped aircraft report
their own identity, position, altitude and velocity roughly once a second. It is a
realistic public stand-in for a sensor feed precisely because it is awkward — fields go
missing, altitude is sometimes the string `"ground"`, positions are stale by a varying
amount, and the same observation is re-served across consecutive polls.

The adapter absorbs all of that. A record without an identity or a position is discarded
rather than passed on degraded. Non-finite and out-of-range values are rejected before
they reach the estimator.

## Temporal normalisation

The feed reports, per aircraft, how old each position already was when the response was
generated (`seen_pos`). **Measured** on the live feed across two consecutive polls: a
median age of 0.31 s, p90 of 3.97 s, a maximum of 48.53 s, and a per-aircraft change in
age between polls ranging −15.76 s to +3.00 s. Of 135 contacts, 17 were re-served
positions that had not moved.

An earlier version of this code stamped every measurement with the tracker's cycle time
and ignored the reported age. That is wrong, and the error is not benign: a constant lag
is invisible to a constant-velocity filter, but the *variation* in age makes a steadily
flying aircraft appear to lurch, and a correctly functioning innovation gate then
rejects good data.

The current design holds one invariant: **a filter's state is valid for exactly one
instant, and that instant advances only when a measurement arrives.**

- A batch carries the instant it arrived; each contact carries the age it reported.
- The step to a new measurement is `(arrival − last_arrival) + last_age − age`. Both
  terms are relative, so no absolute observation time is ever reconstructed and no
  duration is ever subtracted from an `Instant`. There is no underflow path.
- A measurement describing a moment at or before the one the filter already holds is
  counted as *superseded* and dropped. Re-served snapshots land on exactly zero, so
  deduplication falls out of the same comparison.
- The display and the conjunction screen never advance a filter. They take a
  `TrackView`, which extrapolates a copy.

Late measurements are dropped rather than folded back in. Out-of-sequence measurement
handling and retrodiction are **not** implemented; see Future work.

## WGS84 → ENU

Tracking is done in metres in a local tangent frame, not in degrees. Degrees are not a
metric space — at Athens, 0.01° is 1.1 km of northing but 0.87 km of easting — and a
filter run directly on latitude and longitude inherits that distortion.

Geodetic → ECEF → ENU, a rigid rotation with no small-angle approximation. Bowring's
closed form provides the exact inverse.

The inverse matters for a reason worth stating: **ENU "Up" is not altitude.** It is
height above the tangent plane at the observer, and the Earth curves away from that plane
as the square of the range. **Measured**: at 463 km — the edge of the default 250 NM
picture — an aircraft truly at 36,089 ft sits 19,233 ft *below* the Athens tangent
plane. An earlier version displayed that raw value as altitude, an error of ~55,000 ft
at long range which looked plausible at short range. The picture now reports true
altitude above the ellipsoid via the inverse transform.

The tracking and screening frame was deliberately left unchanged. The curvature error is
common-mode between two nearby aircraft: **measured**, a true 305 m vertical separation
at 463 km range computes as 297 m in ENU, an 8 m error against a 305 m threshold. Both
facts are pinned by tests.

## Kalman tracking

Three independent 2-state (position, velocity) filters, one per ENU axis.

A 3-D constant-velocity target with diagonal process and measurement noise has a
block-diagonal covariance — the axes never exchange information — so three scalar
filters are numerically identical to one 6-state filter, with no matrix inversion, no
heap allocation, and a fixed instruction count per update.

Position is initialised from the first plot at the measurement variance. Velocity is
seeded from the reported ground speed and track where present, with a deliberately loose
prior so that two updates overrule the seed. Process noise is the discretised
continuous white-noise-acceleration form.

## Measurement gating

All three axes are evaluated before any is applied — a partially applied update from a
bad plot is worse than a rejected one. A firm track whose innovation exceeds the
configured sigma on any axis rejects the whole plot and counts it.

Coasting to the measurement's own moment happens whether or not the measurement survives
the gate. A rejected plot leaves a coasted track, which is what it should leave.

## Track lifecycle

Tracks are associated by ICAO 24-bit address. That is a property of this sensor, not of
the architecture: ADS-B supplies a unique identity, so no nearest-neighbour association
is required. A radar adapter would need real association, and nothing downstream would
change.

A track becomes *firm* after four agreeing reports; only firm tracks are screened, since
a two-plot track has a velocity that is mostly prior. Tracks are dropped when their last
accepted **observation** — not their last received packet — is older than the configured
timeout. Coasting widens the covariance rather than freezing it, so an aged track
visibly loses confidence in the `±m` column.

## Closest-approach screening

Under a constant-velocity assumption the separation between two tracks is a quadratic in
time, so the closest approach has a closed form and needs no search:

```text
t_cpa = -(dp · dv) / (dv · dv)
```

clamped to `[0, horizon]`, because a closest approach in the past is not a warning. A
cheap reject — can the pair reach the minima within the horizon even closing head-on? —
keeps the O(n²) screen affordable. Both tracks are extrapolated to a common instant
first; comparing two filters at the instants they happen to sit at is how phantom
conflicts are created.

Horizontal and vertical separation are judged separately, against 5 NM and 1000 ft by
default, because that is how airspace is actually divided.

A pair already inside the minima alerts at `t_cpa = 0` even while separating. Loss of
separation now is still loss of separation.

## C ABI

`aether_cpa` exposes the closest-approach primitive to C and C++:

```c
typedef struct { double e, n, u, ve, vn, vu; } aether_state_t;
typedef struct { double t_cpa, horiz_m, vert_m, closing_ms; } aether_cpa_t;

int aether_cpa(const aether_state_t *a, const aether_state_t *b,
               double horizon_s, aether_cpa_t *out);
```

Rules held on that boundary: `#[repr(C)]` on everything that crosses it; every pointer
null-checked before dereference; every input checked for finiteness; `out` untouched on
any error path; no allocation across the boundary; and no panic, since unwinding across
an FFI boundary is undefined behaviour. The math lives in a safe Rust function so it is
unit-tested without an `unsafe` block.

## Input-boundary protection

The ingestion client enforces a **configurable maximum response body of 8 MiB**
(`MAX_BODY_BYTES`). **Measured** justification: three consecutive polls of the default
endpoint returned 91,547 / 91,543 / 91,543 bytes, so the limit is roughly 90× observed
steady state and cannot fire on legitimate traffic even at a much larger radius.

Enforcement is in two places:

1. If `Content-Length` is present and exceeds the limit, the response is rejected before
   any body byte is read. This is advisory only — the header is optional and a hostile
   server can omit or understate it.
2. The body is read chunk by chunk, and the size is checked **before** each copy. A limit
   applied after the copy has already paid for the memory it claims to refuse.

Exceeding the limit drops the response, closing the connection: the remainder is never
read, buffered, or parsed. The JSON parser sits on the success path only, so it cannot
structurally receive more than the limit.

Requests send `Accept-Encoding: identity`. The build cannot automatically decompress —
`reqwest` is configured with `default-features = false` and without gzip, brotli or
deflate — but that is a fact about the manifest, and manifests drift. Asking for identity
makes the assumption explicit at the boundary. **This is not a claim to have addressed
every compression-related attack**; it means the body bound cannot be bypassed by a
response that expands after measurement.

A failed poll produces a batch with no contacts and updated health counters. An empty
contact list is not a claim that the sky is empty; the health line is what distinguishes
the two. Failures are counted separately as HTTP errors, oversized responses, and decode
errors, because an unreachable feed, a flooding feed, and a feed that changed its schema
are three different problems.

**What this protects against:** unbounded memory allocation caused by an oversized or
malformed upstream HTTP response. **What it does not:** anything else. It is not general
network hardening, and Aether has no inbound listener for such hardening to apply to.

---

## Verification

### Tests — 55 passing

These are properties held by construction and checked in CI-able unit tests, not
observations of the live feed. Among them:

- ENU round-trips to geodetic within a millimetre; "Up" is shown to diverge from
  altitude at long range; pair separation is shown to survive that divergence.
- The filter converges on a constant-velocity target and widens its uncertainty when
  coasting.
- Correctly aged measurements pass the gate on a jittered timeline; the same
  measurements with their age discarded are gated. The fixture's jitter is sized to the
  measured feed, and produces a worst innovation of 0.01σ honoured versus 16.1σ
  discarded.
- Viewing a track does not advance its filter.
- Re-served and out-of-order measurements are superseded, not reapplied.
- A gated plot still coasts its track.
- Head-on conflicts are detected; vertically separated and receding traffic is not.
- Tracks last heard at different times are compared at a single instant.
- The body accumulator never exceeds its bound, for chunk sizes from 1 byte to 64 KiB.
- An oversized response is classified as oversized rather than as a transport error, and
  does not terminate the polling loop — verified against a test-only local listener that
  serves 256 KiB with no `Content-Length`.
- The C ABI rejects null pointers and non-finite inputs without dereferencing or
  propagating them.

`cargo clippy --all-targets --all-features -- -D warnings` is clean.

### Controlled A/B — timestamp handling

Pre-fix and corrected binaries, alternating order, 45-second windows separated by
75-second cooldowns, against the same live feed. Alternating order and cooldowns were
necessary because the upstream throttles under load, and throttling otherwise tracks run
order.

| window | build | observations | gated | gated rate |
| ------ | --------- | -----------: | ----: | ---------: |
| 1 | corrected | 2,137 | 3 | 0.14% |
| 2 | pre-fix | 2,752 | 308 | 11.19% |
| 3 | corrected | 2,208 | 81 | 3.67% |
| 4 | pre-fix | 2,755 | 578 | 20.98% |

Pooled: **pre-fix 886 / 5,507 = 16.1% of observations gated; corrected 84 / 4,345 =
1.9%.** Both corrected windows sit below both pre-fix windows, with no overlap.

The pre-fix build also dropped zero tracks in both of its windows, because stamping every
measurement with the arrival time makes every track appear permanently fresh. The
corrected build ages tracks out on the staleness of their evidence and consequently
re-initiates more of them.

### Live baseline — 90 seconds, after all changes

118 tracks · 4,468 observations · 13 gated (0.29%) · 616 superseded · 39 polls ·
2 HTTP errors · 0 oversized · 0 decode errors · mean cycle 0.08 ms · worst cycle 0.24 ms.

Worst cycle is a measurement of this workload on one machine, not a real-time guarantee.
Aether makes no scheduling, priority, or deadline guarantees, and has not been tested
under memory pressure or on constrained hardware.

---

## Security boundary and threat model

**Aether has no inbound production listener.** It is an outbound-only HTTP client. There
is no server, endpoint, route, or socket bind in the shipped binary. `tokio`'s `net`
feature appears only as a dev-dependency, for a test-only listener, and dev-dependencies
are not compiled into the release build.

Consequently there is deliberately **no mTLS, no payload signature scheme, no nonce or
replay infrastructure, and no authentication layer**. Those controls protect an interface
this system does not expose. Adding them would be defending an attack surface that does
not exist.

What is and is not protected:

| concern | status |
| ------- | ------ |
| Upstream server identity, transport integrity | Protected by TLS (rustls) |
| Unbounded memory from an oversized response | Protected by the body limit |
| Non-finite values reaching the estimator | Rejected at the decode boundary and at the C ABI |
| Implausible kinematics | Rejected by the innovation gate |
| Duplicate / out-of-order observations | Detected and counted |
| **Authenticity of an ADS-B observation** | **Not protected. See below.** |
| Physical sensor compromise | Outside the trust boundary |

### ADS-B authenticity

TLS authenticates the upstream server and protects the bytes in transit. It establishes
nothing about whether the physical aircraft state described by an ADS-B message is true.

ADS-B is an unauthenticated broadcast by design. Any transmitter can emit any ICAO
address with any position, and a fabricated observation arrives over a perfectly valid
TLS connection as a perfectly well-formed record. **No cryptography applied at this layer
can fix that**, because the falsehood is introduced before the signed or encrypted
channel begins.

Aether's current defences against fabricated data are behavioural, not cryptographic:
numerical validation, temporal consistency, kinematic gating, and track-state
consistency. These raise the cost of a naive spoof; they do not establish authenticity.
Independent corroboration — multilateration across receivers with known geometry — is the
real defence, and it requires a second sensor this project does not have.

---

## Known limitations

- **Sensor spoofing is unresolved**, as described above. This is the most significant
  limitation in the system.
- **Separation minima are en-route values applied everywhere.** The screen uses 5 NM and
  1000 ft regardless of airspace class, so low-altitude traffic near an aerodrome — where
  much smaller separations are normal and lawful — can raise alerts that a real system
  would suppress. No airspace model is implemented. This is a known false-positive source
  and is not a defect in the closest-approach math.
- **The range column is ENU horizontal distance, not great-circle distance.** The two
  differ by about 0.06% at 463 km (462.3 km against 462.6 km). Not corrected.
- **Altitude is barometric where the feed reports it**, falling back to geometric.
  Barometric altitude is not height above the ellipsoid; the inverse transform corrects
  the reference frame, not the pressure datum.
- **No out-of-sequence measurement handling.** Late observations are dropped and counted,
  not retrodicted.
- **The body limit bounds retained bytes, not the allocator's peak.** `Vec` growth
  doubles, so filling toward the limit can transiently hold both the old and new buffers
  — roughly 12 MiB peak against an 8 MiB bound.
- **The limit bounds one body at a time.** Aggregate memory is bounded in practice
  because there is a single sensor task with one request in flight; that is a property of
  the current architecture, not of the control.
- **A single sensor.** No fusion across sensors of different modalities is implemented,
  despite the architecture being shaped to accept it.
- **Constant-velocity motion model only.** Manoeuvring targets are tracked with elevated
  innovations; no IMM or manoeuvre-adaptive filtering.
- **The tracker has been exercised at ~150 tracks.** Behaviour at thousands is not
  measured, and the O(n²) screen would need spatial partitioning first.

## Future work

- Deterministic capture and replay, so a session can be re-run offline for benchmarking
  and regression testing without depending on the live feed.
- A second sensor modality through the same `Contact` interface — orbital objects
  propagated from public TLE catalogues are the natural next one, and would exercise
  genuine multi-source fusion.
- Airspace-aware separation minima.
- Out-of-sequence measurement handling.
- Spatial partitioning for the conjunction screen.
- Persistence of observations and estimated states for post-hoc analysis.

## Licence

MIT. See [LICENSE](LICENSE).
