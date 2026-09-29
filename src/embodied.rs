//! Embodied-AI boundary types and deterministic goal gating.
//!
//! This module deliberately does not parse natural language and does not drive motors.
//! A language/VLA adapter may ground an utterance into a typed [Goal]. A tracker or
//! perception stack provides a [WorldState]. The deterministic [DecisionGate] then
//! refuses ambiguous, stale or low-confidence targets before producing an [ActionIntent].
//!
//! Keeping these boundaries explicit makes language grounding, tracking and action
//! selection independently testable.

use crate::geo::Enu;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EntityKind {
    Person,
    Robot,
    Vehicle,
    Obstacle,
    Landmark,
    Unknown,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Observation {
    pub source: String,
    pub observed_id: String,
    pub kind: EntityKind,
    pub pos: Enu,
    pub vel: Option<Enu>,
    pub confidence: f64,
    pub age_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EntityState {
    pub track_id: String,
    pub label: Option<String>,
    pub kind: EntityKind,
    pub pos: Enu,
    pub vel: Enu,
    pub confidence: f64,
    pub evidence_age_s: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WorldState {
    pub entities: Vec<EntityState>,
}

impl WorldState {
    /// Resolve by stable track id first, then by a unique human-facing label.
    /// A duplicated label is ambiguous and therefore resolves to None.
    pub fn resolve(&self, target: &str) -> Option<&EntityState> {
        if let Some(entity) = self.entities.iter().find(|e| e.track_id == target) {
            return Some(entity);
        }

        let mut matches = self
            .entities
            .iter()
            .filter(|e| e.label.as_deref() == Some(target));
        let first = matches.next()?;
        matches.next().is_none().then_some(first)
    }
}

/// Auditable output of a language-grounding model. The original utterance is preserved
/// so evaluation can compare the command with the typed goal rather than only the final motion.
#[derive(Clone, Debug, PartialEq)]
pub struct GroundedCommand {
    pub utterance: String,
    pub language: String,
    pub grounding_confidence: f64,
    pub goal: Goal,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Goal {
    Stop,
    Observe {
        target: String,
    },
    Follow {
        target: String,
        distance_m: f64,
    },
    Approach {
        target: String,
        stop_distance_m: f64,
    },
}

/// High-level action only: never a wheel speed, servo position or joint torque.
#[derive(Clone, Debug, PartialEq)]
pub enum ActionIntent {
    Stop,
    Observe {
        track_id: String,
        target_pos: Enu,
    },
    Follow {
        track_id: String,
        target_pos: Enu,
        target_vel: Enu,
        distance_m: f64,
    },
    Approach {
        track_id: String,
        target_pos: Enu,
        stop_distance_m: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RejectReason {
    InvalidGroundingConfidence,
    LowGroundingConfidence,
    MissingOrAmbiguousTarget,
    InvalidEntityConfidence,
    LowEntityConfidence,
    InvalidEvidenceAge,
    StaleEvidence,
    InvalidDistance,
    NonFiniteState,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecisionGate {
    pub min_grounding_confidence: f64,
    pub min_entity_confidence: f64,
    pub max_evidence_age_s: f64,
    pub max_distance_m: f64,
}

impl Default for DecisionGate {
    fn default() -> Self {
        Self {
            min_grounding_confidence: 0.80,
            min_entity_confidence: 0.80,
            max_evidence_age_s: 1.0,
            max_distance_m: 100.0,
        }
    }
}

impl DecisionGate {
    pub fn decide(
        &self,
        command: &GroundedCommand,
        world: &WorldState,
    ) -> Result<ActionIntent, RejectReason> {
        if !unit_interval(command.grounding_confidence) {
            return Err(RejectReason::InvalidGroundingConfidence);
        }
        if command.grounding_confidence < self.min_grounding_confidence {
            return Err(RejectReason::LowGroundingConfidence);
        }

        match &command.goal {
            Goal::Stop => Ok(ActionIntent::Stop),
            Goal::Observe { target } => {
                let entity = self.admit_target(target, world)?;
                Ok(ActionIntent::Observe {
                    track_id: entity.track_id.clone(),
                    target_pos: entity.pos,
                })
            }
            Goal::Follow { target, distance_m } => {
                self.admit_distance(*distance_m)?;
                let entity = self.admit_target(target, world)?;
                Ok(ActionIntent::Follow {
                    track_id: entity.track_id.clone(),
                    target_pos: entity.pos,
                    target_vel: entity.vel,
                    distance_m: *distance_m,
                })
            }
            Goal::Approach {
                target,
                stop_distance_m,
            } => {
                self.admit_distance(*stop_distance_m)?;
                let entity = self.admit_target(target, world)?;
                Ok(ActionIntent::Approach {
                    track_id: entity.track_id.clone(),
                    target_pos: entity.pos,
                    stop_distance_m: *stop_distance_m,
                })
            }
        }
    }

    fn admit_distance(&self, distance_m: f64) -> Result<(), RejectReason> {
        if !distance_m.is_finite() || distance_m < 0.0 || distance_m > self.max_distance_m {
            return Err(RejectReason::InvalidDistance);
        }
        Ok(())
    }

    fn admit_target<'a>(
        &self,
        target: &str,
        world: &'a WorldState,
    ) -> Result<&'a EntityState, RejectReason> {
        let entity = world
            .resolve(target)
            .ok_or(RejectReason::MissingOrAmbiguousTarget)?;

        if !unit_interval(entity.confidence) {
            return Err(RejectReason::InvalidEntityConfidence);
        }
        if entity.confidence < self.min_entity_confidence {
            return Err(RejectReason::LowEntityConfidence);
        }
        if !entity.evidence_age_s.is_finite() || entity.evidence_age_s < 0.0 {
            return Err(RejectReason::InvalidEvidenceAge);
        }
        if entity.evidence_age_s > self.max_evidence_age_s {
            return Err(RejectReason::StaleEvidence);
        }
        if !finite_enu(entity.pos) || !finite_enu(entity.vel) {
            return Err(RejectReason::NonFiniteState);
        }

        Ok(entity)
    }
}

fn unit_interval(x: f64) -> bool {
    x.is_finite() && (0.0..=1.0).contains(&x)
}

fn finite_enu(x: Enu) -> bool {
    x.e.is_finite() && x.n.is_finite() && x.u.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn person(id: &str, label: &str) -> EntityState {
        EntityState {
            track_id: id.to_string(),
            label: Some(label.to_string()),
            kind: EntityKind::Person,
            pos: Enu { e: 4.0, n: 2.0, u: 0.0 },
            vel: Enu { e: 0.5, n: 0.0, u: 0.0 },
            confidence: 0.98,
            evidence_age_s: 0.1,
        }
    }

    fn command(goal: Goal) -> GroundedCommand {
        GroundedCommand {
            utterance: "Ακολούθησε τον Πάνο".to_string(),
            language: "el".to_string(),
            grounding_confidence: 0.97,
            goal,
        }
    }

    #[test]
    fn stop_requires_no_target() {
        let out = DecisionGate::default()
            .decide(&command(Goal::Stop), &WorldState::default())
            .unwrap();
        assert_eq!(out, ActionIntent::Stop);
    }

    #[test]
    fn follow_resolves_unique_label_to_stable_track() {
        let world = WorldState {
            entities: vec![person("person-17", "Panos")],
        };
        let out = DecisionGate::default()
            .decide(
                &command(Goal::Follow {
                    target: "Panos".to_string(),
                    distance_m: 2.0,
                }),
                &world,
            )
            .unwrap();

        match out {
            ActionIntent::Follow { track_id, distance_m, .. } => {
                assert_eq!(track_id, "person-17");
                assert_eq!(distance_m, 2.0);
            }
            _ => panic!("expected follow intent"),
        }
    }

    #[test]
    fn duplicated_human_label_is_not_guessed() {
        let world = WorldState {
            entities: vec![person("person-17", "Panos"), person("person-42", "Panos")],
        };
        let err = DecisionGate::default()
            .decide(
                &command(Goal::Observe {
                    target: "Panos".to_string(),
                }),
                &world,
            )
            .unwrap_err();
        assert_eq!(err, RejectReason::MissingOrAmbiguousTarget);
    }

    #[test]
    fn stale_target_is_rejected() {
        let mut target = person("person-17", "Panos");
        target.evidence_age_s = 3.0;
        let world = WorldState {
            entities: vec![target],
        };
        let err = DecisionGate::default()
            .decide(
                &command(Goal::Follow {
                    target: "Panos".to_string(),
                    distance_m: 2.0,
                }),
                &world,
            )
            .unwrap_err();
        assert_eq!(err, RejectReason::StaleEvidence);
    }

    #[test]
    fn low_language_grounding_confidence_is_rejected_before_motion() {
        let world = WorldState {
            entities: vec![person("person-17", "Panos")],
        };
        let mut c = command(Goal::Follow {
            target: "Panos".to_string(),
            distance_m: 2.0,
        });
        c.grounding_confidence = 0.4;
        let err = DecisionGate::default().decide(&c, &world).unwrap_err();
        assert_eq!(err, RejectReason::LowGroundingConfidence);
    }

    #[test]
    fn non_finite_target_state_is_rejected() {
        let mut target = person("person-17", "Panos");
        target.pos.e = f64::NAN;
        let world = WorldState {
            entities: vec![target],
        };
        let err = DecisionGate::default()
            .decide(
                &command(Goal::Approach {
                    target: "Panos".to_string(),
                    stop_distance_m: 1.0,
                }),
                &world,
            )
            .unwrap_err();
        assert_eq!(err, RejectReason::NonFiniteState);
    }

    #[test]
    fn absurd_or_negative_standoff_is_rejected() {
        let world = WorldState {
            entities: vec![person("person-17", "Panos")],
        };
        for distance_m in [-1.0, f64::NAN, 101.0] {
            let err = DecisionGate::default()
                .decide(
                    &command(Goal::Follow {
                        target: "Panos".to_string(),
                        distance_m,
                    }),
                    &world,
                )
                .unwrap_err();
            assert_eq!(err, RejectReason::InvalidDistance);
        }
    }
}
