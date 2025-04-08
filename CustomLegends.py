import matplotlib.patches as mpatches
import matplotlib.lines as mlines

def LR_handle(lab = 'LR fibers'):
    handle = mpatches.Patch(facecolor='C0', alpha = 0.4, edgecolor='black', label=lab)
    return handle

def HR_handle(lab = 'HR fibers'):
    handle = mpatches.Patch(facecolor='red', alpha = 0.4, edgecolor='black', label=lab)
    return handle

def boundary_handle(lab = None):
    handle = mpatches.Patch(facecolor='None', edgecolor='green', label=lab)
    return handle

def GFA_handle(lab = 'GFAs'):
    handle = mpatches.Patch(facecolor='None', edgecolor='brown', linestyle = '--', label=lab)
    return handle

def safety_margin_handle(lab = None):
    handle = mpatches.Patch(facecolor='None', edgecolor='red', linestyle = '--', label=lab)
    return handle