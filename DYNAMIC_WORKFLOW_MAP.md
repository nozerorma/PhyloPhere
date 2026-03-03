# Dynamic Workflow Map Implementation

This document describes the new dynamic workflow map functionality implemented for PhyloPhere.

## Overview

The workflow map HTML has been completely reimplemented to be **fully dynamic and location-aware**. Instead of relying on workflow parameters to determine completion status, it now scans its environment to detect the presence of canonical PhyloPhere output directories.

## Key Features

### 🔍 **Dynamic Directory Detection**
- Automatically scans parent directory for canonical workflow outputs
- No dependency on original Nextflow execution parameters
- Works regardless of where the HTML file is moved within results structure

### 📱 **JavaScript Refresh Capability**
- "Refresh Status" button updates completion indicators in real-time
- Uses browser fetch API to detect directory existence (where permitted)
- Maintains static generation as primary method with JS as enhancement

### 🚀 **Full Portability**
- Self-contained HTML with embedded CSS and JavaScript
- Can be copied/moved to any location within results directory
- Automatically adapts to new environment on refresh

### 🎯 **Intelligent Stage Detection**
- Uses canonical directory patterns to determine completion status
- Example: `caastools/` directory existence indicates CT workflow completion
- Fallback logic handles partial or incomplete workflow runs

## Implementation Details

### Canonical Directory Mapping

Each workflow stage is detected based on the presence of specific directories:

```groovy
prune        : ['data_exploration/0.Data-pruning']
dataset_rep  : ['data_exploration'] 
ct           : ['caastools', 'discovery']
rer          : ['rerconverge']
fade         : ['selection/fade']
molerate     : ['selection/molerate']
// ... (see WorkflowMap.groovy for complete mapping)
```

### Detection Logic

1. **Primary Detection**: Static scan during HTML generation
   - `scanWorkflowDirectory()` checks for canonical directories
   - `checkStageCompletion()` validates each stage individually
   - Results used to set initial completion status

2. **Secondary Detection**: JavaScript refresh capability
   - Client-side directory detection via fetch API
   - Updates visual indicators without regenerating HTML
   - Handles browser security limitations gracefully

### Integration Points

#### Main Workflow (`main.nf`)
```groovy
workflow.onComplete {
    def outdir = params.outdir ?: "${workflow.projectDir}/Out"
    def ctx = WorkflowMap.buildCtx(outdir, params, workflow)
    // ... generates HTML with dynamic detection
}
```

#### Standalone Generation (`workflows/workflow_map.nf`)
```groovy
GENERATE_WORKFLOW_MAP(trigger) // Uses same dynamic detection
```

#### Manual Generation (`generate_workflow_map.groovy`)
```bash
groovy generate_workflow_map.groovy /path/to/results/
```

## Usage Examples

### Within Nextflow Pipeline
The HTML is automatically generated with dynamic detection:
```bash
nextflow run main.nf --profile docker --outdir ./Results
# Creates ./Results/workflow_map.html with dynamic capabilities
```

### Manual Standalone Generation
```bash
# Generate map for any results directory
groovy generate_workflow_map.groovy /path/to/existing/results/

# Generate in current directory
groovy generate_workflow_map.groovy
```

### Moving/Copying HTML
```bash
# Copy HTML to subdirectory - it will adapt automatically
cp Results/workflow_map.html Results/ct_analysis/
# The copied HTML will detect CT-specific outputs in its new location
```

## Directory Structure Detection

The dynamic map detects completion based on this canonical structure:

```
Results/
├── workflow_map.html           # The dynamic HTML (works anywhere)
├── data_exploration/
│   ├── 0.Data-pruning/        # → prune = COMPLETED
│   └── 2.CT/                  # → contrast = COMPLETED
├── caastools/                 # → ct = COMPLETED  
├── rerconverge/               # → rer = COMPLETED
├── selection/
│   ├── fade/                  # → fade = COMPLETED
│   └── molerate/              # → molerate = COMPLETED
└── HTML_reports/              # Traditional HTML reports (still linked)
    ├── CT_signification.html
    └── RERconverge.html
```

## Browser Compatibility

### Static Detection (Always Works)
- Detection performed during HTML generation
- Status embedded in HTML at creation time
- No browser limitations

### JavaScript Refresh (Browser-Dependent)
- Modern browsers: Full functionality
- Restricted environments: Graceful degradation
- CORS limitations handled automatically

## Testing

### Automated Tests
```bash
# Run comprehensive test suite (requires Groovy/Nextflow)
groovy test_dynamic_workflow_map.groovy
```

### Manual Validation
1. Create test structure with canonical directories
2. Generate workflow map HTML  
3. Open in browser and test refresh functionality
4. Move HTML to different location and verify adaptation

### Example Test Structure
```bash
mkdir -p test_results/{caastools,rerconverge,selection/fade}
groovy generate_workflow_map.groovy test_results/
# Should detect: ct=COMPLETED, rer=COMPLETED, fade=COMPLETED
```

## Migration from Static Version

### Backward Compatibility
- Existing `workflow.onComplete` calls work unchanged
- All publishDir configurations remain valid  
- Traditional HTML report links preserved

### New Capabilities Added
- Dynamic directory scanning replaces parameter-based detection
- JavaScript refresh for real-time updates
- Portable HTML works anywhere in results structure

### Configuration Changes
- No configuration file changes required
- Maintains existing output directory structure
- Compatible with all existing profiles and parameters

## Troubleshooting

### HTML Shows All Stages as "Not Run"
- **Cause**: HTML generated in directory without PhyloPhere outputs
- **Solution**: Move HTML to results directory or regenerate in correct location

### JavaScript Refresh Not Working
- **Cause**: Browser security restrictions or CORS limitations
- **Solution**: Use static detection (always works) or serve via local web server

### Stage Detection Incorrect
- **Cause**: Non-standard directory structure or partial workflow completion
- **Solution**: Check canonical directory mapping in `WorkflowMap.groovy`

### Missing Directories After Workflow
- **Cause**: Workflow failed or didn't run specific stages
- **Solution**: Check Nextflow logs and verify stage completion

## Files Modified

1. **`lib/WorkflowMap.groovy`**: Core implementation
   - Added dynamic directory scanning methods
   - Updated HTML template with JavaScript
   - Modified context building for standalone operation

2. **`workflows/workflow_map.nf`**: Process integration  
   - Updated to use dynamic detection
   - Enhanced documentation

3. **`main.nf`**: Main workflow integration
   - Updated `workflow.onComplete` to use new signature

4. **`generate_workflow_map.groovy`**: Standalone generator (NEW)
   - Command-line tool for manual HTML generation
   - Useful for testing and post-workflow analysis

5. **`test_dynamic_workflow_map.groovy`**: Test suite (NEW)
   - Comprehensive testing of dynamic detection logic
   - Validation of HTML generation and portability