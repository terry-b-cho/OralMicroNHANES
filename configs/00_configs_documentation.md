## Variable Configuration File Format

The variable configuration files (`*_vars.txt`) follow a simple but strict format:

### File Structure
1. **One Variable Per Line**
   - Each line contains exactly one variable name
   - No headers or column names
   - No commas or delimiters
   - No empty lines between variables
   - No comments or annotations

### Example Format
```
VARIABLE1
VARIABLE2
VARIABLE3
```

### Variable Categories:
   - Each file contains variables from a specific category, ordered as follows:
     1. `1_demoWAS_vars.txt`: Demographic variables
     2. `2_oradWAS_vars.txt`: Oral health variables
     3. `3_exWAS_vars.txt`: Dietary and nutritional exposures
     4. `4_pheWAS_vars.txt`: Phenotypic measurements
     5. `5_outWAS_vars.txt`: Health outcomes