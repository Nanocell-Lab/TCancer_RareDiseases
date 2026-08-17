# Rare Oncology 2.0

Módulo que recibe un símbolo génico y detecta **coincidencias proteicas
candidatas** entre:

1. variantes germinales de ClinVar asociadas a una condición con ORPHAcode;
2. mutaciones somáticas observadas en GDC/TCGA.

La implementación reemplaza el `merge` manual de los notebooks por una función
parametrizada, paginada y con deduplicación por caso y mutación. El núcleo no
requiere dependencias externas; los gráficos usan Geneview y Matplotlib como
extras opcionales.

## Instalación

Requiere Python 3.10 o superior.

```bash
cd 2.0
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
from rare_oncology import find_shared_mutations, generate_gene_plots

result = find_shared_mutations(
    "EGFR",
    min_cancer_patients=1,
    program="TCGA",
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
```

Se generan:

- `results/egfr_report.json`: resultado completo y metadatos;
- `results/egfr_shared.csv`: tabla lista para revisión.

Con `--plots`, la subcarpeta `plots/` contiene:

- `*_lollipop.png/.svg`: pista proteica Geneview con variantes, dominios
  UniProt, estructura secundaria y escala de residuos;
- `*_venn.png/.svg`: variantes únicas raras, oncológicas y compartidas;
- `*_recurrence.png/.svg`: pacientes por mutación compartida;
- `*_heatmap.png/.svg`: pacientes por mutación y proyecto TCGA;
- `*_dashboard.png/.svg`: los cuatro paneles en una figura.

El Venn informa cardinalidades exactas, pero usa círculos de tamaño fijo: sus
áreas no son proporcionales al número de variantes.
El umbral `--min-patients` también se aplica al conjunto oncológico de los
gráficos, por lo que el solapamiento visual coincide con las mutaciones
compartidas del reporte.

Los límites de densidad pueden ajustarse con:

```bash
rare-oncology ACVR1 \
  --plots \
  --plot-format png \
  --max-plot-mutations 30 \
  --output-dir results/acvr1
```

Geneview y Matplotlib son dependencias opcionales: el análisis tabular continúa
funcionando sin ellas. Para añadir gráficos a una instalación existente:

```bash
python -m pip install -e '.[plots]'
```

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
- ClinVar y GDC son fuentes vivas; cada reporte conserva la fecha UTC de consulta.

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
