import os
import os.path

from visit import *   # 2

prefix = os.getenv("PREFIX")
pwd   = os.getenv("PWD")
db = "localhost:"+prefix+".xmf"

OpenDatabase(db, 0)
AddPlot("Curve", "operators/Lineout/Fracture", 1, 1)
SetTimeSliderState(TimeSliderGetNStates()-1)
LineoutAtts = LineoutAttributes()
LineoutAtts.point1 = (4.5,0,0)
LineoutAtts.point2 = (4.5,8,0)
SetOperatorOptions(LineoutAtts)
DrawPlots()
xypairs =  GetPlotInformation()["Curve"]
#print xypairs

f = open(os.path.join(pwd,prefix+".v"),"w")
f.write("#x\t v\n")
for i in range(len(xypairs)/2):
    f.write("%e \t %e\n"%(xypairs[2*i],xypairs[2*i+1]))
f.close()
sys.exit()
