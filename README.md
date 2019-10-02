# CCAM

Package to fit a Censored Catch Assessment Model (largely based on the SAM package) and perform a Management Strategy Evaluation.

*under construction*
*ONLY intended for western Atlantic mackerel*

# Installation

devtools::install_github("elisvb/CCAM")

The package requires installation of Rtools and TMB.

# Information

Key differences with the SAM package are:
- Total catches can be censored
- F is separable (to accommodate the previous)
- forecasting can include management procedures
- functions for MSE
- additional plots, tables, etc.

# References

*Model*

- Van Beveren, E., Castonguay, M., Doniol-Valcroze, T., Cadigan, N., Plourde, S., Duplisea, D. (2017). How catch underreporting can bias stock assessment and advice in northwest Atlantic mackerel and a possible resolution using censored catch. Fish. Res. http://dx.doi.org/10.1016/j.fishres.2017.05.015 

*Based on (and gratefully acknowledged)*

- https://github.com/fishfollower/SAM

*Code chunks (censored catch, crl transformation catch)*

- Cadigan, N. 2016a. A state-space stock assessment model for northern cod, including under-reported catches and variable natural mortality rates. Can. J. Fish. Aquat. Sci., 73: 296–308. http://www.nrcresearchpress.com/doi/pdfplus/10.1139/cjfas-2015-0047

- Cadigan, N. 2016b. Updates to a Northern Cod (Gadus morhua) State-Space Integrated Assessment Model. DFO Can. Sci. Advis. Sec. Res. Doc., 2016/022. Centre for Fisheries Ecosystem Research, St. John’s, NL. http://publications.gc.ca/collections/collection_2016/mpo-dfo/Fs70-5-2016-022-eng.pdf

*Information on censored catch*

- Hammond, T. R., and Trenkel, V. M. 2005. Censored catch data in fisheries stock assessment. ICES Journal of Marine Science, 62: 1118–1130. https://academic.oup.com/icesjms/article/62/6/1118/617407/Censored-catch-data-in-fisheries-stock-assessment

- Bousquet, N., Cadigan, N., Duchesne, T., and Rivest, L.-P. 2010. Detecting and correcting underreported catches in fish stock assessment: trial of a new method. Can. J. Fish. Aquat. Sci., 67: 1247–1261.http://www.nrcresearchpress.com/doi/abs/10.1139/F10-051#.WKyhadcrKM8
