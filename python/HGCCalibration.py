import numpy as np

X0="X0"
LI="LI"
rho="rho"
Lead            =  {X0:0.5617,   LI:18.259,   rho:11.35}
Copper          =  {X0:1.435,    LI:15.5141,  rho:8.96}
Brass           =  {X0:1.492,    LI:16.3852,  rho:8.53}
Silicon         =  {X0:9.35,     LI:45.7532,  rho:2.33}
Scintillator    =  {X0:42.1283,  LI:70.13,    rho:1.032}
PCB             =  {X0:17.4097,  LI:48.4238,  rho:1.7}
StainlessSteel  =  {X0:1.735,    LI:16.6,     rho:8.02}
Tungsten        =  {X0:0.35072,  LI:10.3056,  rho:19.3}
Aluminium       =  {X0:8.8751,   LI:38.8623,  rho:2.7}
Foam            =  {X0:441.8,    LI:677.9,    rho:0.0999}

class HGCCalibration:
    ecalMIP = 85./1e6
    hefMIP = 85./1e6
    hebMIP = 1500./1e6
    mips = [ecalMIP,hefMIP,hebMIP]

    def __init__(self,version,useLambda = False,includeShielding = False):
        self.useLambda = useLambda
        if version == "v5":
            #HGCEE
            if includeShielding:
                thermal_shield = self.getWeight(
                        [Aluminium,  Foam,  Aluminium],
                        [2.0,        26.0,  2.0      ])
            else:
                thermal_shield = 0.0

            layer1a = self.getWeight(
                    [Copper,  PCB],
                    [0.5,     1.2])
            layer1b = self.getWeight(
                    [Copper,  Lead],
                    [3.0,     1.0])
            # First section of modules
            mod1a = self.getWeight(
                    [Tungsten,  Copper,  PCB],
                    [1.75,      0.5,     1.2])
            mod1b = self.getWeight(
                    [Copper,  Lead,  Copper],
                    [3.0,     1.0,   3.0   ])
            mod1c = self.getWeight(
                    [PCB,  Copper],
                    [1.2,  0.5])

            # Second section of modules, adding 1.05mm of W and 1.1mm of Lead
            mod2a = mod1a + self.getWeight([Tungsten],[1.05])
            mod2b = mod1b + self.getWeight([Lead],[1.1])
            mod2c = mod1b 
            # Third section of modules, adding 1.4mm of W and 2.3mm of Lead
            mod3a = mod2a + self.getWeight([Tungsten],[1.4])
            mod3b = mod2b + self.getWeight([Lead],[2.3])
            mod3c = mod2b 

            layer30 = layer1a + self.getWeight([Tungsten],[4.2])

            hgcEE = np.zeros(30)
            hgcEE[0] = thermal_shield + layer1a
            hgcEE[1] = layer1b + mod1a

            for l in range(2,9,2): #2-9 
                hgcEE[l] = mod1b #odd
                hgcEE[l+1] = mod1a + mod1c #even
            hgcEE[10] = mod1b

            hgcEE[11] = mod1c + mod2a

            for l in range(12,19,2): #12 - 19 
                hgcEE[l] = mod2b #odd
                hgcEE[l+1] = mod2a + mod2c #even
            hgcEE[20] = mod2b

            hgcEE[21] = mod2c + mod3a

            for l in range(22,27,2): # 22 - 27
                hgcEE[l] = mod3b #odd
                hgcEE[l+1] = mod3a + mod3c #even
            hgcEE[28] = mod3b

            hgcEE[29] = layer30

            #HGCHEF
            hgcHEF = np.zeros(12)
            modhefa = layer1a + self.getWeight(
                    [Brass],
                    [40.])
            modhefb = layer1b
            hgcHEF[0] = modhefa

            if includeShielding:
                hgcHEF[0] += self.getWeight(
                        [StainlessSteel],
                        [15.          ])

            for l in range(1,12):
                hgcHEF[l] = modhefb + modhefa

            hgcHEB = np.ones(12) * self.getWeight(
                    [Brass],
                    [34.5*2+9])
            if includeShielding:
                hgcHEB[0] += self.getWeight( [Copper], [30.])
            
            
            self.weights = np.array([
                hgcEE/hgcEE[0],
                hgcHEF/hgcEE[0],
                hgcHEB/hgcEE[0]])


    def getWeight(self, material_list, thickness_list):
        val = 0
        for m,t in zip(material_list,thickness_list):
            if not self.useLambda:
                val += t/m[X0]
            else:
                val += t/m[LI]
        return val



    def eVtoWeightedMIP(self, edep, hit_type, layer, gen_eta):

        return edep*self.weights[hit_type][layer-1]*np.tanh(gen_eta) / self.mips[hit_type]

