# Listas de genes

## Receptores con actividad quinasa intrínseca

- `receptor_tyrosine_kinases.txt`: 58 receptores tirosina-quinasa humanos de
  las 20 familias IUPHAR/BPS.
- `tgfb_receptors.txt`: receptores serina/treonina tipos I y II de la
  superfamilia TGF-β, más el receptor accesorio TGFBR3.

Estas dos listas cubren las familias humanas de receptores con actividad
**proteína-quinasa intrínseca**. TGFBR3 se conserva por pertenecer al complejo,
aunque no tiene actividad enzimática conocida.

## Receptores vinculados a quinasas

`kinase_linked_cytokine_receptors.txt` contiene 43 subunidades humanas de los
receptores de citocinas tipo I y II catalogados por IUPHAR/BPS en las familias
IL-2, IL-3, IL-6, IL-12, prolactina, interferón e IL-10. No se solapan con las
listas anteriores.

Estos receptores **no son quinasas intrínsecas**: su señalización se asocia
principalmente con quinasas JAK y factores STAT. La lista incluye subunidades
de unión, señalización y receptores señuelo curados dentro de esas familias;
por eso una mutación no implica automáticamente activación de JAK/STAT.

Fuentes de clasificación y nomenclatura, consultadas el 15-08-2026:

- IUPHAR/BPS, familia de receptores de citocinas:
  <https://www.guidetopharmacology.org/GRAC/FamilyDisplayForward?familyId=301>
- IUPHAR/BPS, API REST de familias y genes:
  <https://www.guidetopharmacology.org/webServices.jsp>
- HGNC, símbolos humanos aprobados:
  <https://www.genenames.org/>

Ejemplo de ejecución exhaustiva:

```bash
GENE_LIST="gene_lists/kinase_linked_cytokine_receptors.txt" \
OUTPUT_ROOT="results/kinase_linked_cytokine_receptors" \
DATABASE="results/kinase_linked_cytokine_receptors/rare_oncology.sqlite3" \
./scripts/run_rtk_maximum.sh
```

El runner exhaustivo también solicita por defecto la red humana STRING v12.0,
la guarda en el mismo SQLite y genera `*_string_network.png/.svg`. Puede ajustar
la vecindad y el umbral con `STRING_PARTNERS`, `STRING_SCORE` y
`STRING_NETWORK_TYPE`; las ejecuciones permanecen secuenciales para respetar la
API de STRING.
