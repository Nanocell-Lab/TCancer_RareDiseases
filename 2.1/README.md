# Rare Oncology 2.1

Módulo que recibe un símbolo génico y detecta **coincidencias proteicas
candidatas** entre:

1. variantes germinales de ClinVar asociadas a una condición con ORPHAcode;
2. mutaciones somáticas observadas en GDC/TCGA;
3. evidencia adicional separada por capa: frecuencia poblacional, validez
   gen-enfermedad, observaciones de casos y evidencia oncológica.

La implementación reemplaza el `merge` manual de los notebooks por una función
parametrizada, paginada y con deduplicación por caso y mutación. En 2.1 ClinVar
se amplía con VCV/RCV/SCV mediante EFetch y se añaden conectores opcionales para
AlphaMissense, gnomAD, PanelApp, Orphadata XML, STRING v12.0 y tablas exportadas
de otras fuentes. El núcleo no requiere dependencias externas; los gráficos usan
Geneview y Matplotlib como extras opcionales.

## Instalación

Requiere Python 3.10 o superior.

```bash
cd 2.1
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[plots]'
```

NCBI recomienda identificar la aplicación y limitarse a tres peticiones por
segundo sin API key. El cliente ya aplica ese límite. Es aconsejable definir:

```bash
export NCBI_EMAIL="contacto@institucion.cl"
# Opcional: export NCBI_API_KEY="..."
```

## Uso como función

```python
from rare_oncology import (
    AlphaMissenseEvidenceProvider,
    GnomadEvidenceProvider,
    OrphadataXmlEvidenceProvider,
    PanelAppEvidenceProvider,
    find_shared_mutations,
    generate_gene_plots,
)

result = find_shared_mutations(
    "EGFR",
    min_cancer_patients=1,
    program="TCGA",
    evidence_providers=[
        AlphaMissenseEvidenceProvider(
            database="results/rare_oncology.sqlite3"
        ),
        GnomadEvidenceProvider(max_variants=25),
        PanelAppEvidenceProvider(),
        OrphadataXmlEvidenceProvider("../data/Orphanet/en_product6.xml"),
    ],
)

for mutation in result.shared_mutations:
    print(
        mutation.protein_change,
        mutation.diseases,
        mutation.tcga_projects,
        mutation.cancer_patient_count,
    )

generate_gene_plots(result, "results/plots")
```

La red humana de asociaciones proteicas STRING es otra función complementaria.
Mapea primero el símbolo a un identificador STRING humano, consulta SQLite antes
de la web y conserva una caché exacta por versión, score, tipo de red y número de
vecinos:

```python
from rare_oncology import get_string_network, generate_string_network_plots

string_run = get_string_network(
    "EGFR",
    database="results/rare_oncology.sqlite3",
    required_score=700,
    network_type="functional",
    added_nodes=20,
)
generate_string_network_plots(string_run.network, "results/plots")
```

La búsqueda terapéutica es una función complementaria: usa el gen analizado y,
para la capa CIViC, las mutaciones compartidas ya normalizadas. No cambia el
resultado ClinVar/GDC; ClinicalTrials.gov, Open Targets y la tabla regulatoria
también pueden aportar asociaciones a nivel de gen aunque no haya coincidencias.

```python
from rare_oncology import find_therapy_associations

therapy_result = find_therapy_associations(result)

for association in therapy_result.variant_therapy_evidence:
    print(
        association.protein_change,
        association.therapy_name,
        association.effect,
        association.disease_names,
    )
```

## Historial SQLite y ejecución cache-first

### Functional Git clone

`bootstrap/` is a tracked, sub-10 MiB read-only starter bundle. It provides
stored molecular and therapy examples for EGFR, TGFBR1, ACVR1 and PDGFRB,
available STRING networks, plus normalized FDA/EMA tables. The dashboard loads
the seed automatically and always writes new executions under `results/`.

For direct CLI runs, the tracked regulatory tables can be included with:

```bash
rare-oncology EGFR --therapies \
  --regulatory-file fda=bootstrap/regulatory/fda_biomarkers.csv \
  --regulatory-file ema=bootstrap/regulatory/ema_labels.csv
```

See `bootstrap/README.md` and `bootstrap/manifest.json` for compaction and
provenance details. Rebuild it from complete local snapshots with
`python scripts/build_bootstrap_bundle.py`.

La CLI usa por defecto `OUTPUT_DIR/rare_oncology.sqlite3`. Antes de consultar
ClinVar, GDC o las fuentes terapéuticas busca la ejecución compatible más
reciente. Si no existe o tiene más de 24 horas, consulta las APIs y guarda un
nuevo snapshot; los anteriores permanecen como historial.

```bash
# Primera ejecución: WEB→DB; repeticiones compatibles: HIT
rare-oncology EGFR --therapies --output-dir results/egfr

# Base compartida y vigencia de siete días
rare-oncology EGFR \
  --database data/executions.sqlite3 \
  --cache-max-age-hours 168 \
  --therapies

# Ignorar un snapshot vigente y guardar otra ejecución
rare-oncology EGFR --refresh --therapies

# Prohibir toda consulta web; falla si falta un snapshot vigente
rare-oncology EGFR --offline --therapies

# Comportamiento anterior, sin persistencia
rare-oncology EGFR --no-database
```

La API Python devuelve además si hubo caché:

```python
from rare_oncology import (
    ExecutionDatabase,
    find_shared_mutations_cached,
    find_therapy_associations_cached,
)

run = find_shared_mutations_cached(
    "EGFR",
    database="results/rare_oncology.sqlite3",
    program="TCGA",
)
result = run.result
print(run.cache_hit, run.execution_id)

therapy_run = find_therapy_associations_cached(
    result,
    database="results/rare_oncology.sqlite3",
)

with ExecutionDatabase("results/rare_oncology.sqlite3") as database:
    for execution in database.list_executions(gene="EGFR"):
        print(execution.kind, execution.created_at, execution.cache_hit_count)
```

La clave considera parámetros del análisis, fuentes seleccionadas y huellas de
archivos locales (ruta, tamaño, modificación y SHA-256). No almacena email ni
API key. Con `--offline --plots` tampoco se consulta UniProt: el lollipop usa
una longitud inferida sin dominios ni estructura secundaria.
Los generadores recibidos por la API Python se materializan una sola vez antes
de calcular la clave y ejecutar el análisis. Para proveedores personalizados,
implemente `cache_identity()` y devuelva un diccionario con su endpoint y toda
configuración estable que cambie el resultado; no incluya estado descubierto o
modificado durante `collect()`. Los proveedores incluidos ya cumplen ese
contrato. En sistemas POSIX, la base y sus sidecars se restringen a modo `0600`.
SQLite conserva el JSON completo y no está cifrado: aun así debe protegerse como
los demás resultados de investigación y no usarse como un registro clínico
auditado.

### AlphaMissense con caché por variante

AlphaMissense se consulta opcionalmente mediante Ensembl VEP. Sólo se envían SNV
GRCh38 con alelos A/C/G/T inequívocos y la respuesta se conserva únicamente si
también coincide el gen y el cambio proteico. Los scores y clases aparecen en el
JSON, en columnas propias del CSV de evidencia y en
`alphamissense_scores`, `alphamissense_classes` y
`alphamissense_transcripts` del CSV de coincidencias.

```bash
rare-oncology EGFR \
  --source alphamissense \
  --database results/rare_oncology.sqlite3 \
  --output-dir results/egfr
```

La tabla `alphamissense_variant_cache` usa como clave gen, variante genómica y
versión del dataset. También registra resultados sin predicción; por ello un
análisis posterior del mismo gen consulta SQLite primero y sólo envía a VEP las
variantes nuevas. `--refresh` renueva el análisis general, pero reutiliza esta
caché estable. Con `--no-database` AlphaMissense funciona sin persistencia.
AlphaMissense estima patogenicidad missense, no ganancia o pérdida de función,
y no está validado ni aprobado para uso clínico.

### Inputs estructurales para receptores mutados

La función `generate_variant_structure_inputs` obtiene la secuencia canónica
humana revisada desde UniProt, comprueba que el residuo de referencia coincide y
aplica únicamente sustituciones missense simples. Para cada variante escribe:

- un JSON oficial `alphafold3`, versión 4, para AlphaFold 3;
- un YAML versión 1 para Boltz-2;
- un manifiesto JSON con accession, longitud, SHA-256, archivos y variantes
  omitidas.

AlphaFold 3 no acepta YAML en su implementación local oficial; por eso los dos
formatos no se fusionan en un archivo común.

```python
from rare_oncology import generate_variant_structure_inputs

generated = generate_variant_structure_inputs(
    result,
    "results/egfr/structures",
    variants=["L858R", "T790M"],  # omitir para usar las compartidas
    sequence_database="results/rare_oncology.sqlite3",
    model_seeds=[1, 2, 3],
)
```

Uso por CLI:

```bash
# Todas las sustituciones missense compartidas, máximo 20
rare-oncology EGFR \
  --source alphamissense \
  --structure-inputs \
  --model-seed 1 \
  --model-seed 2 \
  --output-dir results/egfr

# Una variante concreta
rare-oncology TGFBR1 \
  --structure-variant R487W \
  --output-dir results/tgfbr1

# Sin acceso web: usa primero SQLite; si falta, proporcionar FASTA local
rare-oncology EGFR \
  --offline \
  --structure-variant L858R \
  --protein-sequence-file downloads/EGFR.fasta \
  --database results/rare_oncology.sqlite3 \
  --output-dir results/egfr
```

Ejecución posterior de los modelos:

```bash
# Boltz-2 genera el MSA mediante el servidor configurado
boltz predict results/egfr/structures/egfr_l858r_boltz2.yaml --use_msa_server

# AlphaFold 3 ejecuta primero su pipeline de datos/MSA y después inferencia
python run_alphafold.py \
  --json_path results/egfr/structures/egfr_l858r_alphafold3.json \
  --model_dir /ruta/modelos \
  --db_dir /ruta/bases \
  --output_dir results/egfr/alphafold3
```

Con `sequence_database`, la secuencia UniProt se guarda en la tabla
`protein_sequence_cache`; las generaciones posteriores consultan SQLite antes
de usar la web. El generador modela una cadena canónica aislada. No infiere automáticamente
ligandos, dímeros, glicanos, membrana, fosforilaciones ni estequiometría del
receptor; esos componentes deben añadirse de forma explícita según la pregunta
estructural.

El dashboard utiliza el comando enfocado siguiente para generar estos archivos
desde un reporte molecular ya almacenado, sin repetir ClinVar ni GDC:

```bash
PYTHONPATH=src .venv/bin/python -m rare_oncology.structure_cli \
  --report-json results/egfr/egfr_report.json \
  --database results/dashboard_master.sqlite3 \
  --output-dir results/dashboard_runs/EGFR/structures \
  --variant L858R --model-seed 1
```

El resultado conserva tres colecciones independientes:

- `variant_therapy_evidence`: evidencia predictiva aceptada de CIViC para la
  mutación o perfil molecular y, opcionalmente, niveles variante-tumor-terapia
  de OncoKB bajo licencia;
- `regulatory_therapy`: etiquetas/aprobaciones FDA o EMA importadas desde una
  tabla oficial proporcionada por el usuario;
- `trial_and_target_discovery`: ensayos de ClinicalTrials.gov y panorama
  gen-diana-fármaco de Open Targets.

Para restringir el análisis:

```python
result = find_shared_mutations(
    "EGFR",
    projects=["TCGA-LUAD", "TCGA-GBM"],
    min_cancer_patients=2,
    exclude_conditions=["Hereditary cancer-predisposing syndrome"],
    allowed_germline_classifications=["Pathogenic", "Likely pathogenic"],
)
```

## Uso por terminal

```bash
rare-oncology EGFR --output-dir results

# Evidencia abierta adicional y snapshot Orphadata versionado
rare-oncology ACVR1 \
  --source alphamissense \
  --source gnomad \
  --source panelapp \
  --orphadata-file ../data/Orphanet/en_product6.xml \
  --output-dir results/acvr1

# Exportes autorizados o descargados de otras fuentes
rare-oncology EGFR \
  --evidence-file clingen=downloads/erepo.tsv \
  --evidence-file civic=downloads/civic.tsv \
  --evidence-file lovd=downloads/lovd.csv \
  --output-dir results/egfr

# Análisis completo con gráficos PNG y SVG
rare-oncology ACVR1 --plots --output-dir results/acvr1

# Sólo pulmón y glioblastoma; exige recurrencia en al menos dos pacientes
rare-oncology EGFR \
  --project TCGA-LUAD \
  --project TCGA-GBM \
  --min-patients 2 \
  --exclude-condition "Hereditary cancer-predisposing syndrome" \
  --classification ALL \
  --output-dir results

# Reporte terapéutico complementario con las tres capas
rare-oncology EGFR --therapies --output-dir results/egfr

# Seleccionar fuentes abiertas concretas
rare-oncology EGFR \
  --therapy-source civic \
  --therapy-source clinicaltrials \
  --max-therapy-records 100 \
  --output-dir results/egfr

# Añadir OncoKB licenciado sin escribir el token en la línea de comandos
export ONCOKB_TOKEN="su-token-personal"
rare-oncology EGFR \
  --therapy-source oncokb \
  --max-therapy-records 25 \
  --output-dir results/egfr

# Ejecutar toda la familia de receptores TGF-beta (13 genes, tipos I-III)
while IFS= read -r gene; do
  [[ -z "$gene" || "$gene" == \#* ]] && continue
  rare-oncology "$gene" \
    --therapy-source oncokb \
    --output-dir "results/tgfb/${gene}"
done < gene_lists/tgfb_receptors.txt
```

### Otros receptores vinculados a quinasas

Las familias con actividad proteína-quinasa intrínseca ya están cubiertas por
los RTK y los receptores serina/treonina de TGF-β. Como conjunto
complementario, `gene_lists/kinase_linked_cytokine_receptors.txt` contiene 43
subunidades de receptores de citocinas tipo I/II asociados principalmente a
JAK/STAT. No son quinasas intrínsecas y la lista incluye algunas subunidades de
unión o señuelo; consulte `gene_lists/README.md` para el alcance y las fuentes.

```bash
GENE_LIST="gene_lists/kinase_linked_cytokine_receptors.txt" \
OUTPUT_ROOT="results/kinase_linked_cytokine_receptors" \
DATABASE="results/kinase_linked_cytokine_receptors/rare_oncology.sqlite3" \
./scripts/run_rtk_maximum.sh
```

### Ejecución exhaustiva de los 58 receptores tirosina-quinasa

`gene_lists/receptor_tyrosine_kinases.txt` contiene los 58 RTK humanos de las
20 familias IUPHAR. El runner activa registros ClinVar sin filtrar por
clasificación, todos los programas GDC, detalles VCV/RCV/SCV, gnomAD, PanelApp,
CIViC, OncoKB, ClinicalTrials.gov, Open Targets, la red STRING humana y gráficos
PNG/SVG:

```bash
cd 2.1
export NCBI_EMAIL="investigador@institucion.org"
export NCBI_API_KEY="opcional"
export ONCOKB_TOKEN="su-token-personal"

# Opcionales: snapshots locales que amplían enfermedad rara y regulación
export ORPHADATA_FILE="/ruta/en_product6.xml"

# Descarga y normaliza las tablas oficiales FDA/EMA para RTK + TGF-beta.
./scripts/generate_regulatory_tables.py
source downloads/regulatory/regulatory_tables.env

./scripts/run_rtk_maximum.sh
```

Para ignorar snapshots vigentes y consultar nuevamente todas las fuentes:

```bash
REFRESH=1 OUTPUT_ROOT="results/rtk_fresco" ./scripts/run_rtk_maximum.sh
```

Variables operativas: `MAX_LOLLIPOP_RESIDUES` (80), `MAX_LOLLIPOP_LABELS`
(15), `MAX_HEATMAP_MUTATIONS` (30), `MAX_RARE_DISEASE_MUTATIONS` (30),
`MAX_THERAPY_RECORDS` (1000),
`STRING_PARTNERS` (20), `STRING_SCORE` (700), `STRING_NETWORK_TYPE`
(`functional`), `SLEEP_SECONDS` (1), `GENE_LIST`, `OUTPUT_ROOT` y `DATABASE`. Cada gen
se guarda en su propio directorio, los logs quedan en `OUTPUT_ROOT/logs/` y una
base SQLite compartida evita repetir consultas compatibles. `classification
ALL` maximiza recuperación pero mezcla patogénicas, VUS, benignas y conflictivas;
la salida sigue siendo para investigación y requiere revisión de evidencia.

```bash
# Añadir un snapshot regulatorio oficial ya descargado
rare-oncology EGFR \
  --therapy-source civic \
  --regulatory-file fda=downloads/fda_biomarkers.csv \
  --regulatory-file ema=downloads/ema_labels.tsv \
  --output-dir results/egfr
```

Se generan:

- `results/egfr_report.json`: resultado completo, evidencia y estado de fuentes;
- `results/egfr_shared.csv`: tabla lista para revisión.
- `results/egfr_evidence.csv`: registros por fuente y capa, sin colapsar sus
  significados clínicos.
- `results/egfr_therapies.json` y `results/egfr_therapies.csv`: asociaciones
  terapéuticas, sólo cuando se usa `--therapies`, `--therapy-source` o
  `--regulatory-file`.

La capa regulatoria no hace scraping ni infiere aprobación. Cada fila debe
incluir `gene` y `therapy_name`/`drug`; se aceptan, entre otras:

```text
gene, source_record_id, therapy_name, therapy_id, protein_change,
biomarker_text, disease_id, disease_name, therapy_status, effect,
evidence_type, evidence_level, jurisdiction, citation_id, source_url,
source_version, license
```

`therapy_status` debe reflejar el documento original, por ejemplo
`regulatory_approved`, `label_biomarker`, `accelerated_approval`, `withdrawn` u
`orphan_designation`. Una designación huérfana no se convierte en aprobación.

## Generador ejecutable de tablas FDA y EMA

`scripts/generate_regulatory_tables.py` descarga la tabla farmacogenómica de la
FDA y el JSON oficial de páginas de medicamentos de la EMA. Por defecto filtra
los símbolos presentes en `gene_lists/receptor_tyrosine_kinases.txt` y
`gene_lists/tgfb_receptors.txt`, conserva también los snapshots originales y
escribe archivos directamente compatibles con `--regulatory-file`:

```bash
cd 2.1
./scripts/generate_regulatory_tables.py
source downloads/regulatory/regulatory_tables.env
```

Puede usarse con listas propias o con snapshots previamente descargados:

```bash
./scripts/generate_regulatory_tables.py \
  --gene-list gene_lists/receptor_tyrosine_kinases.txt \
  --gene-list gene_lists/tgfb_receptors.txt \
  --genes TP53,KRAS,BRAF \
  --output-dir downloads/regulatory

./scripts/generate_regulatory_tables.py \
  --fda-input snapshots/fda.html \
  --ema-input snapshots/ema_medicines.json \
  --output-dir downloads/regulatory_offline
```

La tabla FDA usa siempre `label_biomarker`: la propia FDA advierte que no todas
las entradas exigen una acción basada en el biomarcador. EMA sólo genera una
fila cuando la indicación terapéutica menciona explícitamente el símbolo o un
alias molecular controlado; conserva por separado estados como autorizado,
condicional, retirado, suspendido o rechazado. Las tablas son evidencia de
investigación y deben contrastarse con la etiqueta o información de producto.

OncoKB se activa sólo con `--therapy-source oncokb`; no forma parte del conjunto
predeterminado de fuentes abiertas. Requiere una licencia válida y
`ONCOKB_TOKEN` (o `--oncokb-token`). El token se envía como encabezado Bearer y
se excluye de la identidad de caché, SQLite y exportes. Los registros conservan
el nivel OncoKB, tumor OncoTree, efecto sensibilidad/resistencia, versión y
PMIDs, pero no redistribuyen las descripciones terapéuticas licenciadas.

## Redes humanas STRING

Para obtener los datos y el gráfico STRING del mismo gen analizado:

```bash
rare-oncology EGFR \
  --string-network \
  --string-score 700 \
  --string-partners 20 \
  --string-network-type functional \
  --plots \
  --database results/rare_oncology.sqlite3 \
  --output-dir results/egfr
```

Para actualizar únicamente STRING (sin repetir ClinVar, GDC ni terapias), use
el ejecutable enfocado que también consulta primero el SQLite:

```bash
PYTHONPATH=src .venv/bin/python -m rare_oncology.string_cli TGFBR1 \
  --database results/dashboard_master.sqlite3 \
  --output-dir results/dashboard_runs/TGFBR1/string \
  --string-score 700 --string-partners 20 --string-network-type functional
```

La ejecución produce:

- `egfr_string_network.json`: red completa y procedencia;
- `egfr_string_nodes.csv`: proteínas/nodos;
- `egfr_string_edges.csv`: scores combinados y canales de evidencia;
- `plots/egfr_string_network.png/.svg`: grafo determinista sin NetworkX.

JSON y CSV conservan todas las aristas. Para evitar un *hairball*, el gráfico
muestra todas las conexiones directas con el gen consultado y sólo las
asociaciones secundarias de mayor score; el pie declara cuántas aristas se
representaron respecto del total exportado.

La especie está fijada a `Homo sapiens` (`9606`). El símbolo se resuelve primero
con `get_string_ids` y después se consulta `network` contra la dirección estable
de STRING v12.0. `--string-score` acepta 0-1000, `--string-partners` 0-100 y
`--string-network-type` distingue asociaciones `functional` y `physical`.
La red se guarda en la tabla `string_network_cache` del mismo SQLite. Una
ejecución posterior compatible devuelve `HIT SQLite`; `--offline` nunca llama a
STRING y falla explícitamente si esa combinación no fue almacenada.

Para una lista de genes, ejecute secuencialmente y deje un segundo entre genes;
no paralelice consultas a la API:

```bash
while IFS= read -r gene; do
  rare-oncology "$gene" \
    --string-network --string-score 700 --string-partners 20 --plots \
    --database results/rare_oncology.sqlite3 \
    --output-dir "results/string/${gene}"
  sleep 1
done < gene_lists/tgfb_receptors.txt
```

Este conector está diseñado para genes analizados individualmente o listas
acotadas. Para el proteoma humano completo debe descargarse el archivo oficial
`9606.protein.links` de STRING y procesarse localmente; la propia documentación
desaconseja miles de consultas a la API.

Una arista STRING funcional puede integrar coexpresión, minería de texto,
bases curadas, experimentos u otras señales. No debe describirse automáticamente
como unión física, causalidad, efecto de mutación ni diana terapéutica.

Con `--plots`, la subcarpeta `plots/` contiene:

- `*_lollipop.png/.svg`: pista proteica Geneview con variantes, dominios
  UniProt, estructura secundaria y escala de residuos;
- `*_venn.png/.svg`: variantes únicas raras, oncológicas y compartidas;
- `*_recurrence.png/.svg`: pacientes por mutación compartida;
- `*_heatmap.png/.svg`: pacientes por mutación y proyecto TCGA;
- `*_rare_disease_heatmap.png/.svg`: asociaciones ClinVar por enfermedad rara
  y cambio proteico, con hasta 30 mutaciones y prioridad para las compartidas
  con cáncer (marcadas con `[C]`);
- `*_dashboard.png/.svg`: los cuatro paneles en una figura.
- `*_string_network.png/.svg`: red STRING humana cuando se incluye
  `--string-network`.

El Venn informa cardinalidades exactas, pero usa círculos de tamaño fijo: sus
áreas no son proporcionales al número de variantes.
El umbral `--min-patients` también se aplica al conjunto oncológico de los
gráficos, por lo que el solapamiento visual coincide con las mutaciones
compartidas del reporte.

El lollipop conserva la coordenada proteica exacta, agrega cambios que caen en
el mismo residuo y reparte el cupo entre variantes compartidas, oncológicas y
raras. Los cambios sinónimos `p.X=` se omiten del gráfico por defecto, aunque
permanecen en JSON/CSV. Los límites de densidad pueden ajustarse por panel:

```bash
rare-oncology ACVR1 \
  --plots \
  --plot-format png \
  --max-lollipop-residues 80 \
  --max-lollipop-labels 15 \
  --max-heatmap-mutations 30 \
  --max-rare-disease-mutations 30 \
  --output-dir results/acvr1
```

El heatmap de enfermedades raras usa el número de registros ClinVar distintos
asociados a cada par enfermedad–cambio proteico. En resultados antiguos que no
conservan ese desglose utiliza presencia/ausencia como compatibilidad. El tope
de 30 es estricto: las mutaciones compartidas con oncología se seleccionan
primero y el cupo restante se completa por abundancia de asociaciones raras.

Use `--include-synonymous-plots` únicamente cuando también quiera representar
cambios proteicos sinónimos. `--max-plot-mutations` se conserva como límite
heredado para comandos anteriores que no indiquen los límites específicos.

Geneview y Matplotlib son dependencias opcionales: el análisis tabular continúa
funcionando sin ellas. Para añadir gráficos a una instalación existente:

```bash
python -m pip install -e '.[plots]'
```

## Análisis global incremental por familias de proteínas

`rare-oncology-global` construye un índice derivado a partir de las ejecuciones
guardadas por 2.1. En la primera ejecución proyecta cada snapshot molecular y
terapéutico a un resumen compacto; después sólo lee filas cuyo `id` sea mayor
que el cursor registrado para cada SQLite. No consulta la web ni recalcula
ClinVar/GDC. La terapia sólo se une cuando su firma molecular coincide
exactamente con el análisis activo del gen, evitando atribuir terapias antiguas
a un conjunto nuevo de variantes.

```bash
cd 2.1
PYTHONPATH=src python -m rare_oncology.global_cli \
  --results-root results \
  --index-database results/global_family_index.sqlite3 \
  --gene-list-dir gene_lists \
  --output-dir results/global_family \
  --plots --plot-format png --plot-format svg
```

La API Python equivalente es:

```python
from rare_oncology import build_global_family_analysis, generate_global_family_plots

report = build_global_family_analysis(
    ["results/rtk/rare_oncology.sqlite3", "results/tgfb/rare_oncology.sqlite3"],
    index_database="results/global_family_index.sqlite3",
    gene_list_dir="gene_lists",
)
generate_global_family_plots(report, "results/global_family/plots", formats=("png", "svg"))
```

Se generan `global_family_analysis.json`, `global_gene_summary.csv`,
`global_family_summary.csv`, `global_family_group_summary.csv` y cuatro figuras:

1. **Matriz familia × capa de evidencia:** porcentaje de solapamientos
   confirmados y porcentaje de genes con evidencia variante-específica,
   biomarcador regulatorio, ensayo o mecanismo/diana.
2. **Paisaje por gen:** variantes raras-oncológicas confirmadas frente a la
   recurrencia máxima en pacientes, coloreado por familia amplia.
3. **Barras apiladas de accionabilidad:** cada gen se asigna únicamente a su
   capa disponible más alta para evitar doble conteo visual.
4. **Heatmap enfermedad rara × subfamilia:** hasta 20 enfermedades abundantes
   y 16 subfamilias, contando asociaciones de mutaciones compartidas con cáncer.

Las familias se toman de las listas versionadas en `gene_lists/` y se subdividen
en grupos biológicamente interpretables (por ejemplo EGFR/ERBB, FGFR, VEGFR,
TGF-β tipo I/II e IL-6). Esta clasificación es analítica, no una inferencia de
causalidad. El campo `actionability_tier` es sólo una vista de ordenamiento: los
conteos crudos de evidencia variante-específica, regulación, ensayos y
mecanismo/diana permanecen separados y no son clínicamente equivalentes.

## Completar consultas ausentes en los SQLite

`rare-oncology-complete` recorre todas las bases 2.1 bajo `results/`, obtiene el
último análisis de cada gen y completa únicamente las capas que no tengan un
estado `ok` o `partial`. No vuelve a consultar ClinVar ni GDC: reconstruye las
identidades GRCh38 almacenadas para añadir AlphaMissense/gnomAD/PanelApp,
fusiona terapias sólo con la misma firma molecular y usa las tablas específicas
de STRING y secuencias UniProt.

Primero se recomienda auditar sin red ni modificaciones:

```bash
cd 2.1
PYTHONPATH=src python -m rare_oncology.completion_cli \
  --results-root results \
  --dry-run \
  --output-dir results/completion
```

Para completar todas las fuentes abiertas predeterminadas:

```bash
PYTHONPATH=src python -m rare_oncology.completion_cli \
  --results-root results \
  --output-dir results/completion
```

Las fuentes predeterminadas son STRING, UniProt, AlphaMissense, gnomAD,
PanelApp, CIViC, ClinicalTrials.gov y Open Targets. También se puede limitar el
backfill, por ejemplo:

```bash
PYTHONPATH=src python -m rare_oncology.completion_cli \
  --results-root results \
  --source string --source uniprot --source alphamissense
```

OncoKB sólo se activa explícitamente y el token no se persiste:

```bash
export ONCOKB_TOKEN='...'
PYTHONPATH=src python -m rare_oncology.completion_cli \
  --results-root results --source oncokb
```

Los archivos regulatorios locales se incorporan con
`--regulatory-file fda=RUTA` o `ema=RUTA`. `--offline` sólo audita/reutiliza
datos presentes, `--retry-errors` vuelve a intentar fuentes con error previo y
`--max-genes N` permite lotes controlados. Cada ejecución genera
`completion_report.json` y `completion_report.csv`; los resultados vacíos o no
disponibles quedan registrados para no repetir consultas inútiles.

## Definiciones y límites

- **Enfermedad rara:** el registro ClinVar tiene una clasificación germinal con
  una referencia explícita a Orphanet. No se usa similitud textual.
- Por defecto sólo se incluyen clasificaciones `Pathogenic`,
  `Likely pathogenic` y `Pathogenic/Likely pathogenic`. Para exploración se puede
  usar `--classification ALL`, pero las VUS o clasificaciones conflictivas no
  deben presentarse como asociaciones causales.
- **Oncología:** ocurrencia somática pública en el programa TCGA de GDC; se puede
  cambiar o eliminar el filtro de programa.
- **Anotación proteica:** el lollipop consulta la entrada humana revisada de
  UniProt seleccionada mediante el símbolo génico exacto. Los dominios se
  representan con las anotaciones `Domain`; las hélices y láminas beta con las
  características `Helix` y `Beta strand` de la proteína canónica.
- La posición dibujada se interpreta sobre la proteína canónica de UniProt. Si
  una variante procede de otra isoforma o transcrito, debe verificarse y
  remapearse antes de interpretar su relación con un dominio. Las posiciones
  superiores a la longitud canónica se omiten y se contabilizan en el título.
- Si no existe una entrada humana revisada o UniProt no está disponible, el
  gráfico marca `longitud inferida (sin UniProt)` y no inventa dominios ni
  estructura secundaria.
- **Común/compartida:** el mismo cambio proteico normalizado aparece en ambas
  fuentes. La salida incluye número de pacientes para diferenciar presencia de
  recurrencia.
- **No es confirmación clínica:** un cambio proteico puede corresponder a más de
  un transcrito o cambio genómico. Antes de interpretar una coincidencia hay que
  confirmar HGVS/transcrito y coordenada/alelos GRCh38.
- En sustituciones simples, el módulo compara también `CanonicalSPDI` de ClinVar
  con el cambio genómico GDC. Esas filas se marcan como
  `confirmed_genomic_and_protein` sólo si todos los registros agregados están
  confirmados. Los agregados parciales se marcan
`mixed_genomic_and_candidate_protein` y preservan cada par ClinVar–SSM en
  `genomic_confirmations`.
- **Capas de evidencia:** una frecuencia gnomAD, una asociación de PanelApp y
  una clasificación ClinVar nunca se fusionan en un único puntaje. El JSON
  conserva `evidence_layer`, procedencia, versión, licencia e identidad genómica.
- **Capas terapéuticas:** CIViC sólo aporta evidencia predictiva curada; Open
  Targets siempre se conserva como `mechanism_target`, incluso cuando informa
  una etapa clínica `APPROVAL`; ClinicalTrials.gov siempre se conserva como
  `clinical_trial` y su alcance como `trial_gene_mention` o
  `trial_variant_mention`, no como biomarcador demostrado. `match_locations`
  indica si la mención apareció en elegibilidad, título, resumen, condiciones o
  intervenciones. Sólo una fuente regulatoria explícita puede producir
  `regulatory_approved`. Si el límite corta resultados del registro, el estado
  de ClinicalTrials.gov queda `partial` en vez de aparentar cobertura completa.
- **Identidad terapéutica de variante:** el conector CIViC enlaza inicialmente
  por cambio proteico y lo marca `candidate_protein_variant`, no como identidad
  genómica exacta. Deben confirmarse transcrito, HGVS y coordenadas/alelos antes
  de trasladar la asociación a una variante individual.
- **ClinVar detallado:** 2.1 consulta el XML VCV después del ESummary y conserva
  RCV, SCV, submitters, revisión, fechas y HGVS. Use `--no-clinvar-details` sólo
  para pruebas rápidas o si EFetch no está disponible. La clasificación y
  revisión de cada RCV se enlazan con su condición mediante `TraitSetID`; los
  SCV/submitters aplanados permanecen como procedencia del VCV y no se atribuyen
  a una condición cuando el XML agregado no demuestra esa relación. Si EFetch
  fue consultado pero no puede resolver un RCV-Orphanet, el estado ClinVar queda
  `partial` y el filtro patogénico no reutiliza la clasificación agregada VCV.
- **gnomAD:** sólo se consulta cuando CanonicalSPDI puede convertirse de manera
  inequívoca a un SNV GRCh38 A/C/G/T. Los indels se marcan como no mapeables y
  no como “ausentes”; requieren normalización VCF antes de habilitarlos. La
  frecuencia poblacional no demuestra patogenicidad. Por seguridad el lote se
  limita a 25 alelos.
- **PanelApp y Orphadata:** describen relaciones gen-enfermedad, paneles y nivel
  de evidencia; no son clasificaciones de una variante individual.
- ClinVar, GDC, gnomAD y PanelApp son fuentes vivas; cada registro conserva la
  fecha UTC de consulta. Para reproducibilidad se recomienda guardar el XML de
  Orphadata y los exportes externos junto con su versión y licencia.

## Formato de tablas externas

`--evidence-file SOURCE=RUTA` acepta CSV, TSV o JSON. Está pensado para
exportes obtenidos legalmente de ClinGen ERepo, DECIPHER, LOVD, CIViC, GenCC,
MSeqDR u otras fuentes con políticas de acceso diferentes. No realiza scraping.
Cada fila debe incluir el gen explícitamente; la API Python sólo permite asumir
un archivo monogénico mediante `assume_gene=True`, dejando esa suposición en la
procedencia.

Las columnas se normalizan sin distinguir mayúsculas y se aceptan alias comunes:

```text
gene, source_record_id, evidence_layer, assembly, chromosome, position,
reference, alternate, transcript, hgvs_c, hgvsp, disease_id, disease_name,
mode_of_inheritance, classification, review_level, evidence_level,
population_frequency, source_version, license
```

Cada fallo opcional queda registrado en `source_statuses` sin perder el análisis
ClinVar/GDC. `--strict-sources` transforma esos fallos en error para pipelines que
exigen cobertura completa.

## Pruebas

Las pruebas usan respuestas controladas y no necesitan red:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Las pruebas de renderizado usan el backend no interactivo `Agg`, verifican la
semántica de cada panel y comprueban que los PNG/SVG sean archivos válidos sin
usar comparaciones frágiles de píxeles.

La prueba EGFR comprueba que una variante ClinVar con ORPHAcode y `L858R` se
detecta cuando GDC reporta la misma alteración, y que se respeta el umbral de
pacientes.

## Prueba real ACVR1 (2026-08-11 UTC)

`rare-oncology ACVR1 --plots --output-dir results/acvr1` se ejecutó contra
ClinVar y GDC Data Release 46.0. Con los filtros patogénicos predeterminados
obtuvo 9 variantes raras, 83 observaciones TCGA y 2 coincidencias:

- `R206H`: 7 pacientes de `TCGA-UCEC`, ORPHA:337 y confirmación genómica;
- `G356D`: 2 pacientes de `TCGA-UCEC`, ORPHA:337 y confirmación genómica.

La ejecución produjo los cinco gráficos en PNG y SVG. Estos resultados son una
validación técnica fechada de fuentes vivas, no una conclusión clínica.

## Smoke test EGFR (2026-08-11 UTC)

La implementación fue ejecutada contra ClinVar y GDC Data Release 46.0:

- con el filtro patogénico predeterminado y excluyendo
  `Hereditary cancer-predisposing syndrome`: 1 variante rara, 425 observaciones
  TCGA y 0 coincidencias;
- en modo exploratorio (`--classification ALL`): 69 variantes raras y 4
  coincidencias con la enfermedad Orphanet 294023: `G917R`, `P644L`, `R831C` y
  `T638M`.

Las cuatro coincidencias exploratorias tenían clasificación `Uncertain
significance` y sólo una observación oncológica cada una. Por tanto demuestran
el funcionamiento técnico del cruce, pero **no evidencia causal o clínica**.

## Fuentes oficiales

- [Acceso programático a ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/)
- [Guía de E-utilities y límites de uso](https://www.ncbi.nlm.nih.gov/books/NBK25497/)
- [GDC Data Analysis API](https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/)
- [Filtros y paginación GDC](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [Geneview mutation tracks](https://github.com/ShujiaHuang/geneview/blob/master/docs/mutation_tracks_guide.md)
- [UniProt REST API](https://www.uniprot.org/help/api_queries)
- [gnomAD v4](https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)
- [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/)
- [Orphadata API y descargas](https://api.orphadata.com/)
- [ClinGen Evidence Repository](https://erepo.clinicalgenome.org/evrepo/)
- [CIViC data releases](https://docs.civicdb.org/en/latest/using/data_releases.html)
- [CIViC GraphQL API](https://docs.civicdb.org/en/latest/api.html)
- [ClinicalTrials.gov API v2](https://clinicaltrials.gov/data-about-studies/learn-about-api)
- [Open Targets GraphQL API](https://platform-docs.opentargets.org/data-access/graphql-api)
- [STRING API v12.0](https://version-12-0.string-db.org/help/api/)
- [Descargas humanas de STRING](https://string-db.org/cgi/download?species_text=Homo+sapiens)
- [FDA pharmacogenomic biomarkers](https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling/)
