# P. furiosus 2,3-BDO Metabolic Network Visualization

Interactive network visualization accompanying the publication on engineered
*Pyrococcus furiosus* strains for 2,3-butanediol production.

🔗 **[View the interactive figure](https://jvjvjvjv.github.io/pfu-bdo-networkvis/)**

**Model repository:** [zhanglab/GEM-iPfu](https://github.com/zhanglab/GEM-iPfu)

## Citation
(placeholder)

## Running the Live App Locally
The version above only supports the default view: BDO-ALS vs Parent-COM, no normalization. To explore other strain comparisons and normalization schemes, launch the full Dash app via Docker:

```bash
docker pull jvjvjvjv/pfu-bdo-networkvis
docker run -p 8050:8050 jvjvjvjv/pfu-bdo-networkvis
```
Then open http://localhost:8050 in your browser.