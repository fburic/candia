from pyopenms import *

RT_LO = 1000
RT_HI = 1050
SWATHS = [623, 639]


mzml_scanfiles = [
    '../experiments/20181120-104433/samples/151204_mzml/151204_YSBN11_mSWATH_29_01-YSBN11.mzML',
    '../experiments/20181120-104433/samples/151204_mzml/151204_YSBN11_mSWATH_29_02-YSBN11.mzML'    
]


for scanfile, i in zip(mzml_scanfiles, range(1, len(mzml_scanfiles) + 1)):
    exp = MSExperiment()
    MzMLFile().load(scanfile, exp)

    spec = []
    for s in exp.getSpectra():
        rt = s.getRT()
        rt = np.around(rt, decimals=4)   
        if RT_LO <= rt <= RT_HI:    

            if s.getMSLevel() == 1:
                spec.append(s)
            else: 
                prec_info = s.getPrecursors()[0]
                prec_window_start = prec_info.getMZ() - prec_info.getIsolationWindowLowerOffset()

                if int(prec_window_start) in SWATHS:
                    spec.append(s)
                    
    exp.setSpectra(spec)
    MzMLFile().store(f'samples/scan{i}.mzML', exp)