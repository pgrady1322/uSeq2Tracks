# uSeq2Tracks Rapid Mode Implementation Summary

## ✅ Completed Features

### 1. Rapid Mode Configuration
- **Added `rapid_mode` parameter** to config.yaml (default: false)
- **Clear documentation** explaining rapid mode vs standard mode
- **Validation logic** to detect rapid mode setting

### 2. Conditional Pipeline Logic
- **Modified rule all** to conditionally skip QC steps when rapid_mode: true
- **Skipped components in rapid mode:**
  - FastQC reports for all samples
  - Adapter trimming summary outputs
  - MultiQC comprehensive report
  - Composite tracks for replicate groups
  - Parameter sweep outputs
  - Genrich additional outputs

### 3. Essential Track Generation
- **Maintained core functionality** in rapid mode:
  - SRA downloads (if needed)
  - Genome indexing
  - Primary pipeline outputs (BigWig tracks, peak calls)
  - UCSC hub files for browser viewing

### 4. Rapid Mode Completion Tracking
- **Added rapid_tracks_complete rule** for completion summary
- **Output directory:** `results/{genome_id}/rapid/`
- **Completion marker:** `rapid_tracks_complete.txt` with summary

### 5. Documentation Updates
- **README.md enhanced** with rapid mode sections
- **Example configurations** for both standard and rapid modes
- **Clear usage instructions** for public dataset processing

## 🎯 Use Cases

### Standard Mode (rapid_mode: false)
- **Research datasets** requiring comprehensive QC
- **Novel samples** with unknown quality
- **Publication-ready analysis** with full documentation
- **Parameter optimization** studies

### Rapid Mode (rapid_mode: true)
- **Public datasets** (ENCODE, TCGA, etc.) with known quality
- **Quick browser track generation** for visualization
- **Streamlined processing** when QC is unnecessary
- **Fast turnaround** for track sharing

## 📁 Output Structure Comparison

### Standard Mode Output
```
results/{genome_id}/
├── qc/
│   ├── fastqc/           # Individual FastQC reports
│   └── multiqc_report.html
├── trimmed/              # Adapter trimming outputs
├── {assay}/              # Assay-specific tracks
├── ucsc/
│   ├── hub.txt
│   └── composite_trackDb.txt
└── genrich_sweep/        # Parameter sweeps
```

### Rapid Mode Output
```
results/{genome_id}/
├── {assay}/              # Essential tracks only
├── ucsc/
│   └── hub.txt           # Basic hub (no composites)
└── rapid/
    └── rapid_tracks_complete.txt
```

## 🚀 Implementation Details

### Configuration Variables
- `RAPID_MODE`: Boolean flag from config
- `RAPID_OUTDIR`: Rapid mode output directory path
- Conditional logic throughout rule all inputs

### Pipeline Behavior
- **When rapid_mode: false**: Full pipeline execution (current behavior)
- **When rapid_mode: true**: Skip QC, generate essential tracks only
- **Output tagging**: All outputs tagged with genome_id regardless of mode

### Backward Compatibility
- **Default behavior unchanged**: rapid_mode defaults to false
- **Existing configs work**: No breaking changes to current setups
- **Progressive enhancement**: Rapid mode is opt-in feature

## ✨ Benefits

1. **Faster processing** for public datasets
2. **Cleaner output** focused on browser tracks
3. **Reduced storage** without QC intermediates
4. **Flexibility** to choose appropriate mode per project
5. **Maintained quality** for research datasets needing full QC

The rapid mode feature successfully provides a streamlined pathway for generating essential genome browser tracks while preserving the comprehensive analysis capabilities of the standard pipeline mode.
