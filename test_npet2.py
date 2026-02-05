from ribctl.lib.npet2.run import run_npet2
import vtk
import pyvista as pv
import warnings
warnings.filterwarnings('ignore')
vtk.vtkObject.GlobalWarningDisplayOff()
pv.set_error_output_file('vtk_errors.log') 

def main():
    ctx = run_npet2("7K00")
    print("run_dir:", ctx.store.run_dir)

if __name__ == "__main__":
    main()
