# CASTRA: Crystallographic Analysis and Structural Transformation Tool

CASTRA is a Python-based GUI software designed for CIF file manipulation, crystal structure analysis, and automated crystallographic computations. It simplifies complex crystallographic tasks such as bond angle calculation, torsion comparison, enrichment ratio analysis, dimer generation, and hydrogen bonding detection through a user-friendly interface.

---

 Features

CIF File Viewer & Editor  
  Load and modify atomic coordinates, symmetry operations, and structure metadata.

Bond Lengths, Angles, and Torsions  
  Automatic computation of geometric parameters with table display.

 Torsion Angle Comparison  
  Compare X-ray derived torsion angles with optimized structures from Gaussian.

 Hydrogen Bond Characterization  
  Identify and visualize potential hydrogen bonds based on van der Waals criteria.

 Dimer Generation from MLC Files  
  Convert monomers to dimers using symmetry operations and export to CIF.

 Enrichment Ratio Analysis  
  Analyze and visualize preferred atomic contacts using data from CrystalExplorer.

AIM Table Generator  
  Analyze QTAIM data from Excel and create bond critical point tables.

 Interactive GUI  
  Built with Tkinter andor PyQt5 for easy navigation.


Installation

1. Clone the repository:

```bash
git clone https://github.com/UdaykumarMani24/CASTRA.git
cd CASTRA


pip install -r requirements.txt

