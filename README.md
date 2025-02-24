Codebase to assess the safety and classification of laser products according to IEC 60825 (specifically, BS EN 60825-1:2014), denoted here as 'the standard'.

tables.py implements various tables from standard, however, only for wavelengths 700-1050 nm at the moment.

spherical_diffuser_default.ipynb is an example how to perform the evaluation (classificaiton, nominal ocular hazard distance (NOHD), skin hazard distance) in line with the standard for pulsed laser light from a spherical diffuser source, for a **default evaluation** (i.e., non-extended source).
  
spherical_diffuser_extended.ipynb is the same as the above but considering the angular subtense of the source, i.e., treating it as an **extended source**.

