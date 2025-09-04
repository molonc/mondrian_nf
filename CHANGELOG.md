# Changelog

## [0.1.9.2] - 2025-09-04
### Changed
- Reordered `BAMMERGECELLS`: now runs earlier, after `CONCATALIGNMETRICS` (previously after `HMMCOPY`).
- Set number of CPUs for `ALIGN` to 1.

## [0.1.9] - 2025-08-01
### Changed
- Split `ALIGN` into separate steps: `ALIGNFASTQSCREEN`, `ALIGNMETRICS`, and `BAMMERGECELLS`.

## [0.1.8]
### Changed
- Increased `maxRetries` from 10 to 20.

## [0.1.7]
### Changed
- CPU and memory adjustments.

## [0.1.6]
### Changed
- CPU and memory adjustments.

## [0.1.5]
### Changed
- CPU and memory adjustments.

## [0.1.4]
### Changed
- CPU and memory adjustments.

## [0.1.3] - 2025-03-17
### Changed
- Removed `numcores` argument from `MUSEQ` call to align with prior changes.

## [0.1.2] - 2025-02-21
### Added
- Support for six optional supplementary references.

## [0.1.1] - 2024-10-31
### Changed
- BCCRC Nextflow fork of Mondrian (from MSKCC).
- Minor column name updates (e.g., `pick_met` vs `cell_call`).
- Performance tuning (cores per process, memory, etc.).
- Environment containers migrated to Azure.

## [0.1.0] - 2024-05-01
### Added
- MSKCC Nextflow version of Mondrian.
