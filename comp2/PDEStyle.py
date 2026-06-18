import matplotlib.colors as mplColors
import matplotlib.pyplot as plt
import numpy as np

## how to use
# %run MyStyle1.py
# or
# import PDEStyle as MS

# MS.ShowMyStyles() # for colors
# ms = MyStyleSet() # object to use like
# ms.air.col # to get the color of 'air'

## toggle the matplotlib colors to these
# ms.SetMPLcolors() 
# ms.SetMPLcolors("off")

class StyleThing:
    def __init__(self, color, linestyle, name = None):
        self.col = color
        self.ls = linestyle
        self.name = name
        return
        
    def printMembers(self):
        style_thing_members = [attr for attr in dir(self) if not attr.startswith('__')]
        print("style object members: \t", style_thing_members)
        return
    
class MyStyleSet:
    def __init__(self, option = "nobles"):
        
        # medium types. these are what populate the other arrays 
        self.air = StyleThing('goldenrod', '-', 'air')
        
        self.lxe = StyleThing('darkviolet', '-', 'lxe')
        self.lxe2 = StyleThing('fuchsia', '-', 'lxe2')

        self.cf4 = StyleThing('deepskyblue', '-', 'cf4')
        # self.cf4_2 = StyleThing('deepskyblue', '-')
        
        self.green = StyleThing('springgreen', '-', 'green')
        
        self.red = StyleThing('red', '-', 'red')
        
        # self.test = StyleThing('#d94ff5', '-')
        
        # line styles
        self.line1 = "-"
        self.line2 = "--"
        self.line3 = "-."
        
        # all colors
        self.names = [x for x in self.__dict__.keys() if isinstance(self.__dict__[x], StyleThing)]
        self.allColsStr = [x.col for x in self.__dict__.values() if isinstance(x, StyleThing)]
        self.allColsObj = mplColors.to_rgba_array(self.allColsStr).reshape(1, -1, 4)
        self.index = -1
    
    def show(self):
        fig, ax = plt.subplots(figsize = (len(self.names)*1.5, 2))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.imshow(self.allColsObj, aspect=0.5)
        for i, (color, name) in enumerate(zip(self.allColsStr, self.names)):
            txt = f"{name} \n {color}"
            ax.text(i, 0, txt, ha='center', va='center')
        plt.show()
        return
    
    ## color cycling ##
    def ResetColorCycle(self, start=0):
        self.index = -1 if start == 0 else start
        return
    def GetNextColor(self, stay = False):
        if stay:
            return self.allColsStr[self.index]
        self.index = (self.index + 1) % len(self.allColsStr)
        colStr = self.allColsStr[self.index]
        return colStr
  
    def SetMPLcolors(self, toggle = "on"):
        if toggle == "off" or toggle == 0:
            # plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.rcParamsDefault['axes.prop_cycle'])
            plt.rcParams.update(plt.rcParamsDefault)
            return
        cycleCols = [getattr(self, name).col for name in self.names]
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=cycleCols)
        return

def ShowMyStyles():
    ms = MyStyleSet()
    print("------------------------------------------------")
    print(ms.names)
    ms.air.printMembers()
    ms.show()
    return

# ------- Other Functions ------- #
def wavelength_to_rgb(wavelength, IR=False):
    '''takes in wavelength in nm and returns an RGB value.'''
    deep_purple = (102/255, 0, 153/255)
    near_infrared = (255/255, 0, 0)
    
    wavelength = float(wavelength)
    if wavelength < 420: ## UV 
        wavelength = 420
    if wavelength > 750: ## Vis
        wavelength = 750
    
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R, G, B = ((-(wavelength - 440) / (440 - 380)) * attenuation), 0, (1.0 * attenuation)
    elif wavelength >= 440 and wavelength <= 490:
        R, G, B = 0.0, ((wavelength - 440) / (490 - 440)), 1
    elif wavelength >= 490 and wavelength <= 510:
        R, G, B = 0.0, 1, (-(wavelength - 510) / (510 - 490))
    elif wavelength >= 510 and wavelength <= 580:
        R, G, B = ((wavelength - 510) / (580 - 510)), 1, 0
    elif wavelength >= 580 and wavelength <= 645:
        R, G, B = 1.0, (-(wavelength - 645) / (645 - 580)), 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R, G, B = (1.0 * attenuation), 0, 0
        
    R = np.clip(R, 0, 1)
    G = np.clip(G, 0, 1)
    B = np.clip(B, 0, 1)
    return R, G, B

def map_value_to_hsv(value, A, B, map = 'hsv'):
    # Get the 'hsv' colormap
    cmap = plt.get_cmap(map)

    # Map the value to the range [0, 1]
    normalized_value = (float(value) - A) / (B - A)

    # Get the color corresponding to the normalized value
    color = cmap(normalized_value)

    return color