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
   - Each file contains variables from a specific category:
     - `6_zimWAS_vars.txt`: Zero-inflated laboratory measurements
     - `3_exWAS_vars.txt`: Dietary and nutritional exposures
     - `4_pheWAS_vars.txt`: Phenotypic measurements
     - `5_outWAS_vars.txt`: Health outcomes
     - `1_demoWAS_vars.txt`: Demographic variables
     - `2_oradWAS_vars.txt`: Oral health variables

### Usage in Parallel Processing
This simple format is particularly useful for parallel processing because:
1. Easy to read and parse
2. No special formatting required
3. Can be directly used in bash scripts with commands 
4. Can be easily combined with Schema Structure Files for parallel analysis: