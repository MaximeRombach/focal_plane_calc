import matplotlib.patches as mpatches
import matplotlib.lines as mlines

def LR_handle(extra_lab = ''):
    handle = mpatches.Patch(facecolor='C0', alpha = 0.4, edgecolor='black', label = '' + extra_lab)
    return handle

def HR_handle(extra_lab = ''):
    handle = mpatches.Patch(facecolor='red', alpha = 0.4, edgecolor='black', label = 'HR ' + extra_lab)
    return handle

def boundary_handle(lab = None):
    handle = mpatches.Patch(facecolor='None', edgecolor='green', label=lab)
    return handle

def GFA_handle(lab = 'GFAs'):
    handle = mpatches.Patch(facecolor='None', edgecolor='brown', linestyle = '--', label=lab)
    return handle

def fiducials_handle(lab = 'Fiducials'):
    handle = mlines.Line2D([0], [0], marker='o', color='w', label = lab,
                          markerfacecolor='red', markersize=6)
    return handle

def safety_margin_handle(lab = None):
    handle = mpatches.Patch(facecolor='None', edgecolor='red', linestyle = '--', label=lab)
    return handle

def final_layout_title(project: str, 
                       vigD:float,
                       nb_robots: int,
                       total_modules: int,
                       total_robots: int,
                       inner_gap, global_gap,
                       out_allowance: float,
                       HR_fibers: int = 0,
                       LR_fibers: int = 0) -> str:

    project_info = "Project: " + project + r" - $\varnothing$" + f"{vigD} mm"
    # project_info = r"$\bf{Project: MUST}$"

    robots_info = f"Total # modules: {total_modules} - Total # robots: {total_robots}"
    modules_info = f" {nb_robots} robots per module"
    if HR_fibers == 0:
        fibers_info = ""
    else:
        fibers_info = f"{LR_fibers} LR fibers - {HR_fibers} HR fibers"

    out_allowance_info = f"Out allowance: {out_allowance * 100} %" # convert to %

    if inner_gap != global_gap:
        figtitle = f"{project_info}\n Semi frameless - {modules_info} \n Inner gap: {inner_gap} mm - Global gap: {global_gap} mm \n {out_allowance_info} \n {robots_info} \n {fibers_info}"
    elif inner_gap == global_gap and global_gap == 0:
        figtitle = f"{project_info}\n Frameless - {modules_info} {robots_info} \n {out_allowance_info} \n {robots_info} \n {fibers_info}"
    else:
        figtitle = f"{project_info}\n Framed - {modules_info} \n Gap: {inner_gap} mm {robots_info} \n {out_allowance_info} \n {robots_info} \n {fibers_info}"

    return figtitle

def module_title(nb_robots: int, module_side_length: float, pitch: float, l_alpha: float, l_beta: float, HR_l_alpha: float, HR_l_beta: float) -> str:
    """Returns the title of the module"""
   
    fiber_arms = f'l_alpha: {l_alpha} mm - l_beta: {l_beta} mm'
    if HR_l_alpha != l_alpha or HR_l_beta != l_beta:
        hr_fiber_arms = f'\n HR: l_alpha: {HR_l_alpha} mm - l_beta: {HR_l_beta} mm'
        fiber_arms += hr_fiber_arms   

    title = f'Workspaces of {nb_robots} robots per module \n Primitive triangle side length: {module_side_length} mm - Pitch: {pitch} mm \n' + fiber_arms

    return title    