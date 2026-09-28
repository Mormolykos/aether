# Aether

A real-time aerospace telemetry engine in Rust: it ingests live ADS-B aircraft
broadcasts, maintains a filtered track picture, and screens every pair of tracks for
loss of separation against ICAO separation minima.

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
- Screens all firm pairs for loss of separation over the look-ahead horizon, in closed
  form, and reports each pair's closest point of approach alongside.
- Exposes the closest-approach primitive over a C ABI for use from C or C++.
- Bounds the memory a single upstream response can cause it to allocate.

It does **not** contain targeting, engagement, weapon, or fire-control functionality of
any kind, and is not intended to.

## Where this problem shows up

Aether is a worked instance of a problem shape rather than a product. The shape is:
a stream of noisy, late, out-of-order position reports arrives from a sensor you do not
own and cannot fix, and something downstream has to make a safety decision on it within
a bounded time budget.

That shape recurs across air traffic control and airspace deconfliction, satellite
conjunction assessment, maritime collision avoidance (AIS in place of ADS-B, hours in
place of minutes), UAS detect-and-avoid, and ground-vehicle sensor fusion. The pieces
carry over more or less directly: geodetic-to-local-frame conversion, one filter per
target, innovation gating against a physically motivated threshold, closed-form closest
point of approach, lifecycle promotion and ageing driven by the age of the evidence
rather than the clock, and a C ABI so the screening primitive can be called from an
existing C or C++ stack without linking Rust into the hot path.

ADS-B was chosen as the sensor because it is public, needs no credentials, and is
genuinely awkward in the ways real sensor feeds are awkward — missing fields, string
sentinels where a number belongs, positions of varying staleness, and the same
observation re-served across polls. Everything reported below is reproducible from this
repository against that live feed.

Aether is not a certified system and makes no airworthiness, operational-safety or
regulatory-compliance claim.

## Running it

```sh
cp .env.example .env
cargo run --release
```

Configuration is read from `.env` in the working directory, or from a path given as the
first argument. It holds the observer position, poll and cycle periods, filter tuning,
separation minima, and the response-body limit.

Every value is checked against a declared domain before the HTTP client or the sensor is
created. The Tokio runtime already exists by then. A value that does not parse as its type
stops the load at that key. A file that parses but cannot be run on is refused with every
offending key named at once, rather than panicking later or quietly switching a check
off. Among the things refused are:

- a zero poll or cycle period;
- a negative, NaN or infinite horizon;
- a NaN, zero or negative gate width;
- negative process noise, or zero measurement variance;
- a zero, negative or non-finite separation minimum;
- a filter tuning, gate, minimum or horizon below the smallest normal double
  (2.2250738585072014e-308), where a double no longer holds the value written;
- an observer position off the globe;
- a feed URL that does not parse, is not `https`, or has no host (`https://`,
  `https:///v2/...`).

Two of the refusals are **operating policy**, not correctness bounds: a body limit above
64 MiB, which is a per-response ceiling and says nothing about total memory; and a poll or
cycle period that is not shorter than the track timeout. Whether the endpoint resolves,
answers or serves the expected JSON is a runtime fact, reported by the sensor-health line
rather than predicted from the text of the URL.

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

The covariance is carried as its Cholesky factor `L`, with `P = L Lᵀ`, never as `P`
itself: square-root filtering after Potter and Bierman, reduced to two states. The
textbook update `P⁺ = (I − K H) P` subtracts nearly equal numbers when a measurement is
far more precise than the prediction. With no process noise, a measurement variance of
1e-20 m² and a 9 ms step, it produced a velocity variance of −2.9e-11 and a position
variance of exactly 0; a variance of 1e-8 m² and a 60 s gap also give the zero. In
factored form, the prediction re-triangularises `[F L | L_Q]` with Givens rotations, and
the update scales the first column of `L` by `√R / √S`. Neither contains that
subtraction, and `L Lᵀ` is symmetric and positive semidefinite for any finite `L`. At
ordinary tunings a test holds it to the textbook filter's numbers, step by step.

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

Horizontal and vertical separation are judged separately, against 5 NM and 1000 ft by
default, because that is how airspace is actually divided. An alert therefore answers one
question: is there any instant in `[0, horizon]` at which the pair is inside **both**
minima at once?

Under a constant-velocity assumption each half has a closed form and needs no search.
Vertical separation is linear in time, so it is inside its minimum on one open interval.
Horizontal separation is judged along the relative track, a straight line. With relative
speed `s = |dv_h|`, the pair sits at `a = dp_h · dv_h / s` along that line and
`d = |dp_h × dv_h| / s` across it, and `d` is the closest the pair ever comes
horizontally. If `d < H`, the pair is inside the minimum while `a + s·t` lies within the
half-chord `w = √((H − d)(H + d))`, one open interval:

```text
( (−w − a) / s ,  (w − a) / s )
```

No position, velocity or minimum is squared on its own. `s` comes from `hypot`, and `a`
and `d` from the unit direction, whose components cannot exceed 1. When the product under
the root is not a normal number, as for minima below about 1e-154, the two square roots
are taken separately. Otherwise the product is used, because `√(x²)` rounds back to
exactly `x`: a pair approaching along its line of centres enters at the exact instant.
An earlier form solved the quadratic
`|dv_h|² t² + 2 (dp_h · dv_h) t + |dp_h|² − H² < 0`. With a configured `H` of 1e-300 m,
`H²` underflowed to 0, and two aircraft in exactly the same place were judged not inside
the minimum.

The pair alerts when the two intervals overlap somewhere in `[0, horizon]`, and the start
of the overlap is the time to loss of separation, shown as `LoS T-`.

Whether they overlap is not decided from the intervals' rounded ends. A window can be
narrower than the spacing of doubles at the instant it happens. Against a 1e-15 m
minimum, two aircraft passing exactly through each other 1 s from now are inside it for
2e-17 s, and both rounded ends of that window are 1.0. The decision is taken instead
from the signs of a few polynomials in the relative state, the minima and the horizon.
Three intervals on a line share a point exactly when every two of them do, and each such
pair is the sign of one polynomial.

- **Every sign is proven.** It is first tried with arithmetic that carries a rigorous
  error bound. When that bound cannot settle it, the sign is computed exactly, as an
  unevaluated sum of doubles.
- **An exact tangent is never classified as a proven breach.** It gives exactly zero,
  which is not a breach. At numerically unresolved scales it may produce an unresolved
  `?` result.
- **Unresolved results can occur even inside the accepted configuration domain,** at
  extreme scales such as a 1e-300 m minimum. Such a pair is reported marked `?` rather
  than dropped. `?` means the arithmetic could not prove breach or non-breach, not that
  the pair is safe.
- **Only the reported `LoS T-` time** is computed from the windows' ends, and it is within
  a rounding of the true instant.

There is no approximate pre-filter in front. Both tracks are extrapolated to a common
instant first; comparing two filters at the instants they happen to sit at is how phantom
conflicts are created.

Alerts are listed soonest loss of separation first. Distance at closest approach only
orders pairs whose loss of separation begins at the same instant.

The 3-D closest approach, `t_cpa = -(dp · dv) / (dv · dv)` clamped to `[0, horizon]`, is
still shown for context, but it decides nothing. It is computed through the unit direction
with no small-speed threshold: a threshold that treated a squared speed below 1e-9 as
parallel reported a pair 1 m apart, closing at 1e-5 m/s, as already at its closest. The C
ABI below calls the same function. Metres of altitude and metres of range
are not interchangeable against minima of 305 m and 9260 m: a pair level at 60 s and 9 km
apart may have its 3-D closest approach at 149 s with 891 m between them vertically. It
loses separation at 57 s, and that is when the screen says it does.

A pair already inside the minima alerts at `LoS T-0` even while separating. Loss of
separation now is still loss of separation.

## C ABI

`aether_cpa` exposes the closest-approach primitive to C and C++:

```c
typedef struct { double e, n, u, ve, vn, vu; } aether_state_t;
typedef struct { double t_cpa, horiz_m, vert_m, closing_ms; } aether_cpa_t;

#define AETHER_OK         0
#define AETHER_ERR_NULL  -1  /* a pointer was null */
#define AETHER_ERR_NAN   -2  /* an input or the horizon was NaN, infinite or negative */
#define AETHER_ERR_RANGE -3  /* a finite input exceeded 1e9 in magnitude, or the
                                result could not be represented */

int aether_cpa(const aether_state_t *a, const aether_state_t *b,
               double horizon_s, aether_cpa_t *out);
```

It is a closest-approach primitive, not a separation verdict. Its `horiz_m` and `vert_m`
are the separations at the instant of least 3-D distance, and a pair can be inside both
minima at a different instant while being outside one of them at that one. A caller that
decides loss of separation from these fields alone repeats exactly that error. Aether's
own screen decides it from the violation intervals above.

Rules held on that boundary:

- `#[repr(C)]` on everything that crosses it.
- Every pointer is null-checked before dereference.
- Every input is checked for finiteness and against the numerical envelope below: each of
  the twelve components within ±1e9. The horizon must be finite and non-negative.
- Every result is checked before it is written, so success carries four finite numbers
  for any input inside that domain. Finite inputs alone did not guarantee that: two states
  of order 1e200 overflowed the squared velocity. With the envelope in place the check is
  defensive; no admitted input reaches it.
- Finite is not the same as accurate. The result is the constant-velocity closest
  approach evaluated in floating point, and, as above, not a separation verdict. It has
  no small-speed threshold, and is tested at relative speeds from 1e-3 m/s down to the
  smallest subnormal.
- `out` is untouched on any error path.
- No allocation crosses the boundary.
- No panic, since unwinding across an FFI boundary is undefined behaviour. The math lives in a safe Rust function so it is
unit-tested without an `unsafe` block.

## Input-boundary protection

The ingestion client enforces a **configurable maximum response body of 8 MiB**
(`MAX_BODY_BYTES`). **Measured** justification: three consecutive polls of the default
endpoint returned 91,547 / 91,543 / 91,543 bytes, so the limit is roughly 90× observed
steady state and cannot fire on legitimate traffic even at a much larger radius. The
configuration refuses a limit above 64 MiB: a bound that can be configured away is not a
bound. 64 MiB is an operating policy, not a derived safe maximum, and it bounds one
response, not the process: parsed contacts, queued batches and tracks are held
separately, and total memory is not bounded or measured here.

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

## Numerical domain

Checking that each input is finite does not keep the estimator finite. A ground speed of
1e308 knots is finite, and so is the velocity it seeds, but multiplying that by four
seconds is not. So values that enter from outside are checked against a declared domain
(`src/domain.rs`) at three boundaries: a contact entering the tracker, together with the
frame it is transformed through; the configuration file; and the C ABI.

- **Contacts are checked at `TrackStore::ingest`**, the only way into the estimator. They
  are not trusted from whichever adapter built them: `Contact` is public, and a library
  caller need not have come through a decoder.
  - Latitude, longitude, altitude and reported age must lie in their domains, or the plot
    is refused and counted as invalid.
  - A reported ground speed, track or vertical rate is converted to SI units and then
    checked. Outside its domain it is treated as not reported, because the position beside
    it is still good.
- **The frame is checked too.** A `Frame` records when it is built whether its observer
  position lies in the domain, and ingest refuses every contact transformed through one
  that does not. A finite frame at an altitude of 1e200 m used to place tracks at
  −1e200 m.
- **A filter update is computed on copies and kept only if every state and covariance
  term is finite and no velocity component exceeds 1e9 m/s.** Otherwise the plot is
  refused as invalid and the track is left exactly as it was. That holds for the invalid
  outcome only: a gated plot deliberately coasts its track, and gated and superseded plots
  are counted.
- **The gate is written so that only a finite sigma inside it passes**, and residuals
  that are not numbers are refused before the gate on every track, young or firm. Under
  the old form, "reject if sigma exceeds the gate", a NaN compared false and was
  accepted.
- **The covariance is valid by construction**, symmetric and positive semidefinite,
  because it is carried as a Cholesky factor (see Kalman tracking).

Each quantity has its own bound and its own reason, listed in `src/domain.rs`. Two kinds of
number appear there:

- **1e9 is a policy envelope**, applied to altitude, age, speed and climb components,
  filter tuning, gate width, separation minima, horizon and the C ABI components, in each
  one's SI unit. It is set far above the regional air picture Aether is built for and is
  the same number in every unit only for simplicity. It is not a physical limit: a
  measurement variance above 1e9 m², a sigma above 31.6 km, is not impossible, and this
  demonstrator declines it. Velocity is bounded per axis, so a speed can exceed 1e9 m/s.
- **The smallest normal double, about 2.2e-308, is a derived floor** for measurement
  variance, gate width, separation minima, and non-zero process noise and horizon. Below
  it a double holds fewer than 53 significant bits, so the value used is not the value
  written, and a variance that small has no representable covariance even when its
  factor does. Above it, the filter update and the horizontal interval are formed
  without the squares and subtractions that failed at small values, and are tested down
  to the floor. How small a minimum or gate is sensible is operating policy, and is not
  decided here.

What this does not do:

- It does not bound intermediate arithmetic. A tuning of q = 4 with r = 1e-300, a
  1e-150 s step and a 10 km residual give a candidate velocity near 1e154 m/s. The
  commit check refuses it; it is still computed.
- Transformed positions are not held to the envelope: an observer and a contact each
  inside it can be about twice the envelope apart.
- It is not a geometric operating domain: the envelope bounds magnitudes, not how far
  from the observer the tangent-plane frame stays meaningful.
- The Kalman filter type itself does not validate its inputs; the store is the boundary.

---

## Verification

### Tests — 96 passing

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
- A pair inside both minima at any instant of the horizon alerts, including when its 3-D
  closest approach lies outside the breach, and a pair inside both minima now alerts at
  equal velocity.
- Across 4,000 generated geometries — aimed at the cylinder's rim, at equal and barely
  differing velocities, and at a zero horizon — the screen agrees with an oracle that
  only samples positions every 0.05 s: every sampled breach is alerted, and at every
  alert's onset the pair is on or within a micrometre of both minima, with no earlier
  sampled instant inside both. A breach shorter than the sampling step is invisible to
  that oracle.
- An earlier loss of separation always ranks ahead of a later one, whatever the
  distances.
- Tracks last heard at different times are compared at a single instant.
- The body accumulator never exceeds its bound, for chunk sizes from 1 byte to 64 KiB.
- An oversized response is classified as oversized rather than as a transport error, and
  does not terminate the polling loop — verified against a test-only local listener that
  serves 256 KiB with no `Content-Length`.
- The C ABI rejects null pointers and non-finite inputs without dereferencing or
  propagating them.
- The C ABI either succeeds with four finite numbers or fails without writing a byte of
  the caller's output. This is checked with a bit-pattern sentinel across magnitudes up to
  1e308 in every field and horizons up to `f64::MAX`, including the finite 1e200 states
  whose squared velocity overflowed.
- A pair 1 m apart closing at 1e-5 m/s has its closest approach at the end of a
  100,000 s horizon, not now; the same holds at speeds down to 1e-12 m/s, and one
  geometry scaled from 1e-100 down to the smallest subnormal gives the same answer.
- The covariance stays valid in the case that broke the textbook update: no process
  noise, a measurement variance of 1e-20 m², a 9 ms step. So it does across measurement
  variances from 1e-8 m² to the smallest double, steps from a microsecond to a minute,
  200 stationary updates each, and 60 generated runs of 2,000 noisy steps, where the
  factor is checked after every step and the matrix wherever its entries are
  representable. At ordinary tunings the filter reproduces the textbook filter's state
  and covariance to 1e-9 relative.
- After that very precise report, identical stationary reports are still accepted, and
  the position uncertainty never reads zero.
- Two tracks in exactly the same place breach any positive minimum, down to the smallest
  subnormal, standing still, moving together or moving apart. The same holds through a
  configuration that accepts such a minimum.
- A breach narrower than the spacing of doubles is found: 100 m apart, closing at
  100 m/s, against a 1e-15 m minimum. At that scale, a tangent is not a breach, and one
  ulp outside is not a breach either. One ulp inside is.
- Two windows both narrower than that spacing are told apart exactly. When their centres
  are closer together than the spacing, they overlap or not by their widths alone.
- A tangent in a direction whose unit vector no double holds exactly is still not a
  breach.
- Every decision is the same with the pair's order reversed.
- The proven signs agree with exact integer arithmetic on 20,000 cases, tangencies
  included, at scales from 2^-900 to 2^900.
- In 40 generated runs, no track went non-finite. The runs mix ordinary plots with NaN,
  infinities, `f64::MAX`, 1e200 and subnormals in every field. They come in through the
  public ingest path, which bypasses the decoder, and every track is viewed up to a
  million seconds on. That is evidence over this corpus, not a proof over every
  sequence.
- A plot whose arithmetic would produce a velocity beyond the envelope is refused. A
  frame outside the domain lets no contact in, whether it is non-finite or finite (an
  altitude of 1e200 m, a latitude of 1e10); the track is left untouched.
- Every numeric configuration key is swept through NaN, infinities, zero, negatives and
  extremes. Each value is either refused by name or runs the tracker and the screen
  without a panic or a non-finite number. Subnormal tunings, minima and horizons are
  refused; the smallest normal double and 1e-300 are accepted. Feed URLs that are plain
  `http://`, `ftp://`, missing the colon, missing a host (`https://`, `https:///v2/...`)
  or containing a space are refused by name.

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

### Live baseline — 90 seconds, measured before the September 2026 numerical changes

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
| Upstream server identity, transport integrity | Protected by TLS (rustls). The configuration refuses a feed URL that does not parse as `https` with a host, and the client is built `https_only`, which in reqwest refuses any non-TLS hop, redirects included. A caller that builds the sensor directly can still give it `http` |
| Unbounded memory from an oversized response | Protected by the body limit |
| Values outside the numerical domain | Refused at the tracker's ingest, whoever the caller, together with the frame they are transformed through, and at the C ABI; a filter update is kept only if every result is finite and every velocity component within 1e9 m/s |
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
  and is not a defect in the separation math.
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
