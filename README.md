# KCWI Sky Wizard
**Version: 0.1.2024-12-19**

---

## About this tool

The KCWI Sky Wizard is designed to refine:
- **Sky subtraction**
- **Flux calibration**
- **Telluric correction**

for KCWI/KCRM data. Users must first reduce the data using the [KCWI Data Reduction Pipeline](https://kcwi-drp.readthedocs.io/en/latest/).

### Key Features:
- **Sky Subtraction**: Utilizes the [`ZAP`](https://zap.readthedocs.io/en/latest/) package with a PCA approach.
- **Telluric Correction**: Leverages the `tellfit` function of the [`PypeIt`](https://pypeit.readthedocs.io/en/release/telluric.html) package to derive a best-fit telluric model of the standard star.

---

## Installation

1. **Create a new conda environment:**

   ```bash
   conda create -n kskywizard python=3.11
   conda activate kskywizard
   ```

2. **Clone this repository and install the required packages:**

   ```bash
   git clone https://github.com/zhuyunz/KSkyWizard.git
   cd KSkyWizard
   conda env update --file environment.yml
   ```

3. **Install PypeIt:**

   ```bash
   pip install pypeit
   ```

4. **Install `zap_for_kcwi` in a separate directory:**

   ```bash
   git clone https://github.com/jasonpeng17/zap_for_kcwi.git
   cd zap_for_kcwi
   pip install .
   ```

5. **Return to the KSkyWizard directory, and run the setup code:**
  At the moment, we encourage you to install the KSkyWizard under development mode, as it would be easier to track all the updates.

   ```bash
   cd /path/to/KSkyWizard
   pip install .
   ```

## Download the Mauna Kea telluric data

To perform telluric correction, PypeIt uses preconstructed PCA models. Follow these steps to download and configure the data:

1. **Download the telluric data:**

   ```bash
   pypeit_install_telluric TelFit_MaunaKea_3100_26100_R20000.fits
   ```

2. **Configure the telluric data location.**

   - The telluric data is typically stored in:
   
      `~/.pypeit/cache/download/url/`.
      
   - Each subdirectory contains:
      - A `contents` file (the actual data file)
      - A `url` file (an ASCII file linking to the download URL)
      - Take a note of the full path of the `contents` file, e.g., `~/.pypeit/cache/download/url/5f17ecc1fcc921d6ec01e18d931ec2f8/contents`

   - Locate the installed version of the `processing.cfg` file. Note that this is typically different from where you 
   stored cloned version of this repository. 
      
      ```bash
      kskywizard --config_path
      ```

   This will return the full path of the `procecssing.cfg` file. Modify the `telgridfile` variable in the file, and point it to the content path from above.
   
   Example:

   ```plaintext
   telgridfile = ~/.pypeit/cache/download/url/5f17ecc1fcc921d6ec01e18d931ec2f8/contents
   ```

   **Important**: **DO NOT** include quotation marks around the file path.


## Prerequisites

Run the KCWI Data Reduction Pipeline:
   
   - Process the data using [`KCWI DRP`](https://kcwi-drp.readthedocs.io/en/latest/), preferably following the [`KCWIKit` instruction](https://github.com/yuguangchen1/KcwiKit?tab=readme-ov-file#instuctions).

   - (Optional) Skip the sky subtraction since this tool handles it, by turning on the `-k` flag for the science frames.
   
      Eample:

      ```bash
      reduce -r -f kr*fits -k -g
      ```

      After the science frames, perform sky subtraction for the standard stars to obtain the preliminary sensitivity curves.
      
      For example:
   
      ```bash
      reduce -r -f [Standard Star Frame].fits
      ```

      Alternative, create a `std.list` that contains all the standard star frames and run
      
      ```bash
      reduce -r -l std.list
      ```



## Usage

### 1. Start KSkyWizard

   You can run this tool either in the directory containing the reduced data or from another directory, specifying the input directory later in the GUI.
   
   ```bash
   conda activate kskywizard
   kskywizard
   ```

   A GUI window will appear.
   <img width="1200" alt="image" src="https://github.com/user-attachments/assets/c297825f-c714-4b33-96a2-98eb615cedda">

### 2. Set the input and output directory

- **Input Directory**: Browse and select the `redux` directory where all KCWI DRP outputs are stored.
- **Output Directory**: Browse and select an output directory. **Tip**: Use a separate working directory to avoid overwriting the DRP output.

### 3. Flux calibration and telluric correction.

This process is divided into two steps:
1. **Refining the Sensitivity Curve**: Adjust the sensitivity curve to remove the effects of telluric absorption and sky lines. This step involves interactive region selection and refinement.
2. **Fitting the Telluric Features**: Use the `PypeIt` PCA model to fit the telluric features and finalize the calibration.

#### Step 1: Load and Refine the Sensitivity Curve
1. **Open the DRP Sensitivity File**: 
   - Browse and select the `*_invsens.fits` file from the DRP outputs using the `Open DRP invsens` button.
   
   ![Select DRP Invsens File](https://github.com/user-attachments/assets/9a3afb1e-26dc-43ae-8576-b629e50c32b1)

2. **Interactive Refinement**:
   - To refine the sensitivity curve, click on the plotting canvas to activate the interaction mode.
   - Use the following keyboard shortcuts to adjust regions:
     - **`e`**: Exclude a region from the fit.
     - **`i`**: Include a region in the fit.
     Adjusting the regions will change the automatically generated spline knots, which can be further adjusted by
     - **`a`**: Add a spline knot.
     - **`d`**: Delete a spline knot.
   
   - **Re-fit the Sensitivity Function**:
     - Press **`f`** to re-fit the sensitivity function iteratively until the model (cyan line) aligns well with the data.
     
     ![Interactive Refinement](https://github.com/zhuyunz/KcwiKit/assets/33030986/a0125d01-9fe0-4bc0-807f-d857b89b0463)

   - **Use the DRP’s Default Curve**:
     - Press **`b`** to revert to the DRP's original sensitivity curve if you prefer the DRP version.

#### Step 2: Fit the Telluric Features
1. **Derive the Telluric Model**:
   - Press **`t`** to generate the telluric correction model. This process may take some time, and the GUI will become unresponsive during the calculation. Wait until the output appears.

   ![Telluric Model Output](https://github.com/zhuyunz/KcwiKit/assets/33030986/14e5ab89-2ac4-42b4-9c01-01be296a35b4)

2. **Save the Results**:
   - Once satisfied with the sensitivity curve and telluric fit, press the `Save updated invsens` button to save the refined data.
   These updated curves can be reloaded to the tool using the `Open updated invsens` button for future usage.

---

#### Tips:
- **What to Include and Exclude During Fitting?**  
  - For the **sensitivity curve fitting** (light green region):
    - **Include**: Only the stellar continuum.  
    - **Exclude**: 
      - Residual cosmic rays.
      - Stellar absorption lines (e.g., H-alpha).  
      - Sky emission lines.  
      - Telluric features.  
  - For the **telluric fitting** (dark green region):
    - **Include**: Known telluric features to ensure accurate modeling. For a summary of these features, refer to [Buton et al., 2013](https://www.aanda.org/articles/aa/full_html/2013/01/aa19834-12/aa19834-12.html).

- **Iterative Refinement**:  
  - All fitting can be adjusted and rerun until reaching a satisfactory result. 

### 4. Sky subtraction

   a.  Switch to the science panel and enter the science image number. For instance, enter 85 and press the `return ` button if the frame is kr230923_00085. If you want to use the off-field sky frame, please check the `Use Off-field Sky Frame No.`, also enter the sky image frame number and press return. You can also update the science image number by pressing `Previous` and `Next` button.
   
   b. If this is the first time that you reduce this frame, please press `Load Raw Cube`, `Save Cropped Cube` and `Load Cropped Cube`. Otherwise, simply press `Load Cropped Cube` to load the data. This step is simply to crop the datacube outside of the good wavelength region to avoid running into edge problem with ZAP.

   c. If you plan to use the in-field sky, you need to mask out the science target. You can check the `_wlimg.fits` in the output directory,  mask the source using DS9, and save the DS9 region as kr230923_00085.reg for frame kr230923_00085. After that, enter the mask frame number, press enter, and press `Update ZAP mask` button. You can update ZAP mask multiple times whenever you press the button

   d. Now it's time to run ZAP. 
   - If you want to use the default ZAP setup, please uncheck `Use multiple skyseg in ZAP` and `Additional Sky Seg near Halpha`, and press `Run ZAP`. In this case, it adopts one sky segment and fit the overall variation across the entire wavelength. By default, it adopts `cfwidth=300` to remove the science signal. For more information related to the `cfwidth`, please check the github [page](https://github.com/jasonpeng17/zap_for_kcwi/blob/master/doc/index.rst) of the modified version and the original ZAP paper by [Soto et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3210S/abstract). To proceed, please type 'q' in the text widget like the example below and press enter to start ZAP.
     <img width="1160" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/e320095f-625b-48dc-8bc8-301f4f439044">
     
     Alternatively, you can update the cfwidth by changing the text directly (the last value in the tuple), and type 'u' to update. The screenshot below shows an example to update the cfwidth to 100.
     <img width="1160" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/48e7def7-e149-4c6f-8e71-4ae3c8ea8db4">
     
     The text widget will then print out the updated cfwidth. If you are happy with it, type 'q' in the text wdiget to proceed. If you want to reset the default setting, you can type 'r' instead.
   - We find that using multiple sky segment sometimes can reduce the sky residual. This approach was adopted in the original ZAP by [Soto et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3210S/abstract) when it first came out, but the updated version switched to one sky segment. For how the sky segment is chosen, please refer to [Soto et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3210S/abstract) for more details. If you want to use this option, please check `Use multiple skyseg in ZAP`. Besides, we find that if a very strong emission line happens to fall in a sky segment with relatively sparse sky lines (e.g., Halpha emission), a large `cfwidth` would aritificially introduce more noise in the sky-subtracted spectrum. By checking `Additionally Sky Seg near Halpha` and setting the redshift, the code would generate an additional segment in the region of Halpha with `cfwidth=5`.
     <img width="1149" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/fd3cb55d-6935-47e2-88e0-267ee55428ed">

     If you are happy with this segment, type 'q' in the text widget to proceed. Otherwise, you can update the segment in the text widget directly to whatever regions you want to use, and type 'u' to update the change. Every tuple indicates (the start of the segment, the end of the segment, cfwdith). You can either delete a tuple or add a new one, or simply changing some of the values in it. The screenshot below shows an example when you want to combine segment (7200, 7700, 300) and (7700, 8265, 300) into a new segment with cfwidth = 50.
     <img width="1157" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/51f3fc57-47d7-425c-8dd0-768468c5718b">

     Now the text widget would print out the updated sky segment. Type 'q' if you wanna proceed, or 'r' if you want to reset the sky segment.
     <img width="1149" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/116ee095-48f3-46a1-8c3b-237b24494e8e">

   - If you uncheck or check the box mistakenly, you can simply update the check box and press `Run ZAP` button again, which will reactivate the text widget.
   - It will take a while for ZAP to run, please be patient. The output will be in the terminal window.

### 5. Examine the output
   Now you can examine the sky-subtracted spectrum to see if that looks good to you. By default, it will extract the spectrum from a 10x10 box centered on the center of the FoV (i.e., a box with lower left at (x1, y1) = (13,44) and upper right at (x2, y2) = (23, 54). You can update the region used to extract the spectrum simply by updating the DS9 pixel indices in the box below the plotting canvas as `x1, y1, x2, y2` (lower left + upper right). 
   <img width="1178" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/527c7005-6f49-4b5e-be97-4bfa37c6728c">

   The pre-ZAP spec is scaled to have the same median as the sky-subtracted spectrum for better display purposes.

   If you set the redshift, it will also indicate where the emission lines are expected to be. You can use the navigation bar to examine the plot, and go back to the full plot with right click. Besides, you can also examine the plot with the following options (click the canvas first to activate the keyboard interaction):
   
   - Press 'n' to turn off the pre-ZAP spec.
   - Press 's' to turn on the pre-ZAP spec so that you would know where the skylines are.
   - Press 'o' to move the spectrum to the observed frame. It is useful when you want to update the sky segment (which requires in the observed frame) and re-run ZAP.
   - Press 'r' to move the spectrum to the rest frame.

   Go to the next frame if you are happy with the spectrum; otherwise, press `Run ZAP` button and run the sky subtraction again with the updated sky segment. It's always a good idea to check the output `_zap_icubes.fits`! For the frame using the in-field sky, you can also check the white-lighted image of the clean datacube `_zapclean_wlimg.fits` to make sure the source is properly masked. You can update the DS9 region mask and re-run the ZAP again!


## Contact

If you encounter any problems or find bugs, please email [Zhuyun Zhuang](mailto:zzhuang@astro.caltech.edu).
