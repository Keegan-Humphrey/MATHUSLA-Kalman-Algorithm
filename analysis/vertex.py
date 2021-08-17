import ROOT as root
import numpy as np
import detector


##H FILE
tracking_file_name = "../tracker_files/jan23/h/stat0.root"
tracking_file = root.TFile.Open(tracking_file_name)
tree = tracking_file.Get("integral_tree")


tracking_file_namew = "../tracker_files/jan23/w/stat0.root"
tracking_filew = root.TFile.Open(tracking_file_namew)
treew = tracking_filew.Get("integral_tree")


c = 29.97
wTot = 0.
wTotal = 0.
WinLayer = 0.
WinLayer_1sigma = 0.
WinLayer_2sigma = 0.
WinLayer_3sigma = 0.

hTot = 0.
hTotal = 0.
HinLayer = 0.
HinLayer_1sigma = 0.
HinLayer_2sigma = 0.
HinLayer_3sigma = 0.

det = detector.Detector()

for k in range(int(tree.GetEntries())):
	tree.GetEvent(k)

	if (tree.NumVertices < 1):
		continue

	hTot += 1.
	
	
	y = tree.Vertex_y[0]

	if y > 8900. or y < 6110.:
		continue

	hTotal += 1.

	ey = tree.Vertex_ErrorY[0]

	if det.inLayer(y) >= 0:
		HinLayer += 1.

	if det.inLayer_w_Error(y, ey) >= 0:
		HinLayer_1sigma += 1.

	if det.inLayer_w_Error(y, 2.*ey) >= 0:
		HinLayer_2sigma += 1.

	if det.inLayer_w_Error(y, 3.*ey) >= 0:
		HinLayer_3sigma += 1.





for k in range(int(treew.GetEntries())):
	treew.GetEvent(k)

	if (treew.NumVertices < 1):
		continue

	wTot += 1.

	y = treew.Vertex_y[0]

	if y > 8900. or y < 6110.:
		continue

	wTotal += 1.
	
	ey = treew.Vertex_ErrorY[0]

	if det.inLayer(y) >= 0:
		WinLayer += 1.

	if det.inLayer_w_Error(y, ey) >= 0:
		WinLayer_1sigma += 1.

	if det.inLayer_w_Error(y, 2.*ey) >= 0:
		WinLayer_2sigma += 1.

	if det.inLayer_w_Error(y, 3.*ey) >= 0:
		WinLayer_3sigma += 1.



print("w:")
print(wTotal/wTot)
print([WinLayer/wTotal, WinLayer_1sigma/wTotal, WinLayer_2sigma/wTotal, WinLayer_3sigma/wTotal])

print("h:")
print(hTotal/hTot)
print([HinLayer/hTotal, HinLayer_1sigma/hTotal, HinLayer_2sigma/hTotal, HinLayer_3sigma/hTotal])