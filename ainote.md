# ainote.md: shared progress notes between workstreams

Each workstream keeps its own section below (Done / In progress / Next / Blocked-on).
Please also add a 2–3 sentence **plain-language note** whenever you finish something, so the beginner guide (GUIDE.md) can be updated to match what was actually built. Record any bug the verification caught (needed for the AI-use statement).

## Guide (GUIDE.md)

- **Done:** First draft of GUIDE.md from the plan: what is being built, each piece in plain language (heat PDE, finite elements, time stepping, Newton/nonlinearity, KS aggregation, verification ladder, adjoint, SLSQP, KKT multipliers, results.json + website), glossary, rubric map, AI-use statement skeleton.
- **In progress:** Nothing.
- **Next:** Update each "Status" box and fill the *TBD* numbers as the other workstreams report here (solver/ladder results, gradient checks, optimum and multipliers, viewer). Fill the AI-use table and "bug the verification caught".
- **Blocked on:** Other workstreams' plain-language notes and results.
- **Needs from others:** (1) Final mesh sizes and braking scenario numbers. (2) Ladder table with errors. (3) One bug the ladder or gradient check caught. (4) Active constraints and multiplier values at the optimum. (5) With vs without nonlinearity comparison.
