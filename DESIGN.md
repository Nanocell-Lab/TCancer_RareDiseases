# Design

## Source of truth
- Status: Active
- Last refreshed: 2026-08-17
- Primary product surfaces: Streamlit global landing page and per-gene explorer.
- Evidence reviewed: `dashboard/app.py`, `dashboard/.streamlit/config.toml`, `dashboard/README.md`, dashboard figures/transforms/catalog/runner modules and tests.

## Brand
- Personality: scientific, calm, evidence-forward and explicit about uncertainty.
- Trust signals: source provenance, immutable execution identifiers, evidence-layer separation and visible research-use caveats.
- Avoid: clinical recommendation language, decorative genomics imagery, alarmist colors and dashboard density without hierarchy.

## Product goals
- Goals: summarize the current SQLite corpus before gene-level inspection; reveal family-level convergence; make missing coverage visible; preserve a direct path to exact gene evidence.
- Non-goals: patient-specific treatment advice, causal inference from STRING, or recomputing molecular evidence in the presentation layer.
- Success signals: users can identify high-signal families and genes, understand evidence coverage, and open an exact stored snapshot without a web request.

## Personas and jobs
- Primary personas: translational researchers, molecular oncologists and rare-disease analysts.
- User jobs: assess corpus coverage, compare receptor families, find shared rare-disease/oncology variants, inspect actionability layers and drill into one gene.
- Key contexts of use: exploratory review of an evolving local SQLite corpus and research-meeting preparation.

## Information architecture
- Primary navigation: `Global overview` and `Gene explorer`.
- Core routes/screens: corpus landing page; family/actionability analysis; per-gene molecular, therapy, STRING and structure tabs.
- Content hierarchy: scope and freshness → global metrics → family comparisons → gene ranking → exact snapshot provenance.

## Design principles
- Evidence before decoration: every aggregate must resolve to stored genes and snapshots.
- Progressive disclosure: global overview first, exact records and provenance on drill-down.
- Comparable but not conflated: variant evidence, regulatory status, trials and mechanism evidence remain separate.
- Tradeoffs: favor readable bounded charts and deterministic top-N selections over rendering every point simultaneously.

## Visual language
- Color: slate foundations, orange primary accent, blue rare-disease evidence, purple oncology evidence and green confirmed/actionable evidence.
- Typography: Streamlit sans-serif with compact uppercase eyebrow labels and strong numeric hierarchy.
- Spacing/layout rhythm: 8 px base rhythm, 14–18 px card radius and generous separation between analytical sections.
- Shape/radius/elevation: low-elevation white analytical cards; stronger gradient only for the page hero.
- Motion: native Streamlit transitions only; no decorative animation.
- Imagery/iconography: semantic symbols such as stars for shared variants; no generated imagery.

## Components
- Existing components to reuse: hero, metric cards, Plotly figures, Streamlit tabs, dataframes, download buttons and provenance expanders.
- New/changed components: global scope banner, family filters, family landscape, actionability matrix, convergence heatmaps and ranked gene table.
- Variants and states: populated, filtered, empty corpus, incompatible snapshot and partial-coverage states.
- Token/component ownership: theme tokens remain in `.streamlit/config.toml`; component-specific CSS remains in `app.py` until a second Streamlit page requires extraction.

## Accessibility
- Target standard: WCAG 2.1 AA where Streamlit/Plotly permit.
- Keyboard/focus behavior: all filters and drill-down controls use native Streamlit widgets.
- Contrast/readability: never encode evidence only by color; include labels, symbols and hover text.
- Screen-reader semantics: titles and captions explain chart purpose; tables provide the complete values represented by plots.
- Reduced motion and sensory considerations: no required animation or flashing content.

## Responsive behavior
- Supported breakpoints/devices: desktop-first, usable on tablet; mobile receives stacked metrics and columns through Streamlit's native layout.
- Layout adaptations: charts use responsive width and bounded heights; dense matrices expose filters and top-N caps.
- Touch/hover differences: essential values remain in tables and labels rather than hover alone.

## Interaction states
- Loading: spinner only while invoking the external 2.1 producer.
- Empty: explain that no compatible analysis snapshot exists and point to Gene explorer.
- Error: show source/database error without hiding already readable snapshots.
- Success: show corpus scope and most recent stored execution.
- Disabled: disable gene-run actions when no valid symbol is selected.
- Offline/slow network: render SQLite content immediately; web retrieval remains an explicit gene-level action.

## Content voice
- Tone: concise, technical and cautious.
- Terminology: `shared mutation`, `confirmed genomic and protein match`, `therapy evidence layer`, `stored snapshot` and `protein family`.
- Microcopy rules: distinguish observations from interpretations; use `evidence` rather than `recommendation`; identify top-N truncation.

## Implementation constraints
- Framework/styling system: Python, Streamlit and Plotly; no additional frontend framework.
- Design-token constraints: preserve the current orange/slate theme and existing category colors.
- Performance constraints: read cached SQLite snapshots; do not call web sources from global rendering; derive each global summary once per selected snapshot set.
- Compatibility constraints: dashboard remains an independent package and invokes 2.1 only through its CLI.
- Test/screenshot expectations: pure aggregate and figure tests first; Streamlit startup smoke test; no pixel snapshots.

## Open questions
- [ ] Confirm whether future family catalogs will remain curated text files or become stored 2.1 content.
- [ ] Decide whether Plotly point selection should eventually navigate directly to a gene without an explicit selector.
