# Layered world memory

Embodied agents should not treat every visible thing as equally temporary.

Aether's prototype now separates world state into three layers:

- **Structural** — walls, fixed geometry, columns and other anchors that normally do not move.
- **Persistent movable** — tables, chairs, cabinets and other objects that are stable for long periods but can move.
- **Dynamic** — people, animals, carried objects and short-lived scene contents.

This state lives outside the language/VLA model. It is not a Transformer KV cache.

## Why not use KV cache as world memory?

KV cache accelerates attention over tokens that have already been processed. It is an
inference optimization, not a reliable semantic database. When observations change,
selectively invalidating or editing only the relevant attention state is model-dependent
and can be unsafe.

Instead, keep the physical world as explicit typed state. The model or planner receives
only the local/relevant slice plus the changes since its last revision.

## Update model

Each world object carries a monotonic revision. Consumers can ask for:

- the static map once;
- persistent movable objects when they change;
- live dynamic tracks continuously;
- only objects changed since revision N.

This gives the behavior Panos described: a wall need not be reintroduced every frame,
while a newly observed glass or a moving person can update immediately.
