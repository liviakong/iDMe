from coffea.hist import Cat, Bin

############################
##### Categorical Axes #####
############################
# General Purpose
sample = Cat("sample","Sample Name")
cut = Cat("cut","Cut Applied")
# For electrons
ele_type = Cat("ele_type","Electron Type")
ele_set = Cat("set","Electron Collections Used in Matching")
scheme = Cat("scheme","Matching Scheme")
# For vertices
vtx_sign = Cat("sign","Vertex sign (q1*q2)")
vtx_type = Cat("type","Vertex Type")


############################
###### Numerical Axes ######
############################
# Electrons
ele_pt = Bin("pt","$p_{T}$ [GeV]",100,0,50)
ele_eta = Bin("eta","$\\eta$",60,-3,3)
ele_phi = Bin("phi","$\\phi$",64,-3.2,3.2)
ele_trkHits = Bin("numTrkHits","Number of Tracker Hits",30,0,30)
ele_pixHits = Bin("numPixHits","Number of Pixel Hits",10,0,10)
ele_stripHits = Bin("numStripHits","Number of Strip Hits",25,0,25)
ele_chi2 = Bin("chi2","Track $\\chi^2/df$",200,0,100)
ele_trkIso = Bin("trkIso","Tracker Iso",100,0,100)
ele_trkRelIso = Bin("relIso","Tracker Relative Iso",50,0,5)
ele_prob = Bin("prob","Electron Track $\\chi^2$ Probability",100,0,1)
ele_angRes = Bin("angRes","Angular Resolution $\\sqrt{\\sigma_\\eta^2 + \\sigma_\\phi^2}$", 100, 0, 0.1)
ele_dxy = Bin("dxy","Electron Track $d_{xy}$",50,0,5)
ele_dxySignif = Bin("dxy_signif","Electron Track $d_{xy}/\\sigma_{d_{xy}}$",150,0,150)

# Vertices
vxy = Bin("vxy","$v_{xy}$ [cm]",2000,0,50)
vxy_coarse = Bin("vxy","$v_{xy}$ [cm]",100,0,50)
vxy_zoom = Bin("vxy","$v_{xy}$ [cm]",200,0,20)
vxy_zoomzoom = Bin("vxy","$v_{xy}$ [cm]",500,0,5)
vz = Bin("vz","$v_{z}$ [cm]",2000,0,50)
mass = Bin("mass","$m_{e^+e^-}$ [GeV]",100,0,50)
energy = Bin("energy","Dielectron Energy [GeV]",200,0,100)
dR = Bin("dR","$\\Delta R_{e^+e^-}$",80,0,4)
vtx_vxyErr = Bin("sigma_vxy","$v_{xy}$ Error",100,0,10)
vtx_chi2 = Bin("chi2","Vertex Fit $\\chi^2/df$",150,0,30)
vtx_vxySignif = Bin("vxy_signif","Vertex $v_{xy}$ Significance",100,0,100)
vtx_prob = Bin("prob","Vertex $\\chi^2$ Probability",100,0,1)
vtx_pt = Bin("pt","Selected Vertex $p_T$",100,0,100)

# Gen-matching
etype = Bin("Etype","Electron Match Type",4,0,4)
ptype = Bin("Ptype","Positron Match Type",4,0,4)
nEmatch = Bin("nEmatch","Number of Gen Electron Matches",5,0,5)
nPmatch = Bin("nPmatch","Number of Gen Positron Matches",5,0,5)
matchClass = Bin("matchClass","Match Class",4,0,4)

# Vertex Gen-Matching
vtx_matchType = Bin("type","Vertex Type",7,0,7)

# Jet / MET
jet_met_dPhi = Bin("jet_met_dphi","$\\Delta \phi(\\mathrm{jet},E_T^\mathrm{miss})$",32,0,3.2)