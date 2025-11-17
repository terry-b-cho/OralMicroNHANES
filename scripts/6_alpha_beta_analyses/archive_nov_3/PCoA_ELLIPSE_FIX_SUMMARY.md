# PCoA Ellipse and Centroid Fix Summary

## Problem Identified
PCoA plots were failing to display:
1. 95% confidence ellipses - not rendered at all
2. Centroids - + symbols not colored properly (all black)

## Root Cause Analysis

### Comparison: Working Code vs Broken Code

**Working Code (5_age_analyses_full.R):**
```r
p_pca_integrated <- ggplot(pca_df_full, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  geom_point(alpha = 0.65, size = 0.2) +
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.3, level = 0.90) +
  geom_point(data = pca_centroids, aes(x = PC1_center, y = PC2_center),
             size = 1, shape = 3, stroke = 0.5, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors)
```

**Broken Code (Original):**
```r
p_pcoa <- ggplot(pcoa_data, aes(..., group = !!sym(cat_var))) +  # ❌ Explicit group breaks stat_ellipse
  stat_ellipse(geom = "polygon", ...) +
  stat_ellipse(geom = "path", ...) +  # ❌ Redundant second call
  geom_point(data = centroids, ..., inherit.aes = FALSE, color = "black")  # ❌ Prevents color inheritance
```

## Key Issues Found

1. **Explicit `group` parameter**: The `group = !!sym(cat_var)` in the main `aes()` was interfering with `stat_ellipse`'s automatic grouping mechanism
2. **Dual ellipse calls**: Two separate `stat_ellipse` calls (polygon + path) were conflicting instead of complementing
3. **Centroid color inheritance broken**: `inherit.aes = FALSE` and `color = "black"` prevented centroids from inheriting the proper group colors

## Fixes Applied

### Fix 1: Remove Explicit Group Parameter
- **Before**: `aes(x = PC1, y = PC2, color = !!sym(cat_var), fill = !!sym(cat_var), group = !!sym(cat_var))`
- **After**: `aes(x = PC1, y = PC2, color = !!sym(cat_var), fill = !!sym(cat_var))`
- **Why**: ggplot2 automatically groups by the `color` and `fill` aesthetics. Explicit `group` was interfering with `stat_ellipse`

### Fix 2: Single Ellipse Call
- **Before**: Two calls - `stat_ellipse(geom = "polygon", ...)` + `stat_ellipse(geom = "path", ...)`
- **After**: Single call - `stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.4, level = 0.95, type = "t", segments = 200)`
- **Why**: The working code uses a single ellipse call. The `linetype = "dashed"` parameter creates the dashed outline within the polygon geom itself

### Fix 3: Centroid Color Inheritance
- **Before**: `geom_point(data = centroids, aes(...), inherit.aes = FALSE, color = "black")`
- **After**: `geom_point(data = centroids, aes(x = PC1_center, y = PC2_center, color = !!sym(cat_var)), show.legend = FALSE)`
- **Why**: Centroids need to inherit the color aesthetic to match their group. The centroids data frame already has the categorical variable column from `group_by()`, so it will properly map to colors

## Final Fixed Code

```r
# Calculate centroids - IMPORTANT: Keep the grouping column so centroids can inherit colors
centroids <- pcoa_data %>%
  group_by(!!sym(cat_var)) %>%
  summarise(PC1_center = mean(PC1, na.rm=TRUE), PC2_center = mean(PC2, na.rm=TRUE), .groups = "drop")

# Create PCoA plot with ellipses and centroids
p_pcoa <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = !!sym(cat_var), fill = !!sym(cat_var))) +
  geom_point(alpha = 0.15, size = 0.2) +
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.4, level = 0.95, type = "t", segments = 200) +
  geom_point(data = centroids, aes(x = PC1_center, y = PC2_center, color = !!sym(cat_var)),
             size = 1, shape = 3, stroke = 0.5, show.legend = FALSE) +
  scale_color_manual(values = plot_colors) +
  scale_fill_manual(values = plot_colors, guide = "none") +
  ...
```

## Expected Results

After these fixes:
1. ✅ 95% confidence ellipses will be displayed with dashed outlines around each group
2. ✅ Centroids (+ symbols) will be colored according to their group using `plot_colors`
3. ✅ All three PCoA plots (Braycurtis, Unwunifrac, Wunifrac) will render correctly

## Testing

Running full analysis pipeline:
- Mode 1: Normal full mode
- Mode 2: no_tick_and_bracket_asterisk
- Mode 3: no_tick_and_bracket

All modes use HALF plot style with full dataset (N=9,349 samples).

