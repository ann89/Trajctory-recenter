import math
import MDAnalysis
import numpy as np
#import matplotlib.pyplot as plt


top = 'pbc-mol-4-6us.pdb'
traj = 'pbc-mol-4-6us.xtc' 
#read the vectors
u = MDAnalysis.Universe(top,traj)

def getCMx(xpos, Lx):
   
    PI= math.pi
    
    radius_i = Lx/ (2. * PI)
    
    x_cm_i = np.mean(radius_i * np.cos(xpos / Lx * 2. * PI), axis=0)
    z_cm_i = np.mean(radius_i * np.sin(xpos / Lx * 2. * PI), axis=0)

    return (np.arctan2(-z_cm_i, -x_cm_i) + PI) / (2. * PI) * Lx

def set_pos_back_to_box(res_idxs, CMx, Lx, pos_x_all_new):
    pos_x_all_new = pos_x_all_new.copy()
    for i in range(len(res_idxs)):
        idx = res_idxs[i]
        pos_res = u.select_atoms('resnum %i'%idx).center_of_geometry()
        if pos_res[0] -  CMx + (0.5 * Lx) >= Lx:
            #print('True')
            idx_atoms = u.select_atoms('resnum %i'%idx).ix
            pos_x_all_new[idx_atoms] = pos_x_all_new[idx_atoms] - Lx 
        elif pos_res[0] -  CMx + (0.5 * Lx) < 0:
            idx_atoms = u.select_atoms('resnum %i'%idx).ix
            pos_x_all_new[idx_atoms] = pos_x_all_new[idx_atoms] + Lx
            #print('second')
    return pos_x_all_new


trj = "DSPC_recetered.xtc"  
with MDAnalysis.Writer(trj, multiframe=True, bonds=None, n_atoms=u.atoms.n_atoms) as XTC:
    

    for ts in u.trajectory:
        pos_x_all = u.select_atoms('all').positions[:,0]
        #pos_DAPC_x = u.select_atoms('resname DAPC').positions[:,0]
        pos_DSPC_x = u.select_atoms('resname DSPC').positions[:,0]
    
        Lx = u.trajectory.ts.triclinic_dimensions[0][0]

        #CMx = getCMx(pos_DAPC_x[:],Lx)
        CMx = getCMx(pos_DSPC_x[:],Lx)
        pos_x_all_new = pos_x_all -  CMx + (0.5 * Lx)

        DAPC_idxs = u.select_atoms('resname DAPC and name P').resnums
        pos_x_all_new = set_pos_back_to_box(DAPC_idxs, CMx, Lx, pos_x_all_new)

                
        DSPC_idxs = u.select_atoms('resname DSPC and name P').resnums
        pos_x_all_new = set_pos_back_to_box(DSPC_idxs, CMx, Lx, pos_x_all_new)
            

                
        CHL_idxs = u.select_atoms('resname CHL and name O2').resnums
        pos_x_all_new = set_pos_back_to_box(CHL_idxs, CMx, Lx, pos_x_all_new)
            

    
       # plt.figure()
       # plt.scatter(pos_x_all_new, u.select_atoms('all').positions[:,1]) 
       # plt.scatter(pos_x_all_new[16490:16547], u.select_atoms('all').positions[16490:16547,1]) 
       # plt.show()
        
        #trj = "new_trajectory.xtc"  
        #with MDAnalysis.Writer(trj, multiframe=True, bonds=None, n_atoms=u.atoms.n_atoms) as XTC:
            
    
       
        foo1 = np.asarray([pos_x_all_new, u.select_atoms('all').positions[:,1], 
                                       u.select_atoms('all').positions[:,2]])
        foo2=foo1.transpose()
    
        u.atoms.positions = foo2
        
        XTC.write(u.atoms)
                
            
#     #u.select_atoms('all').positions[:,0] = u.select_atoms('all').positions[:,0] - CMx
#     pos_x_all_new = pos_x_all -  CMx + 0.5 * Lx
#     pos_x_DAPC_new = pos_DAPC_x - CMx + 0.5 * Lx
#     pos_x_others_new = pos_x_others - CMx + 0.5 * Lx
#     pos_x_all_new[pos_x_all_new>Lx] = pos_x_all_new[pos_x_all_new>Lx] - Lx
#     
#     print(CMx, Lx)
#     plt.figure()
#     plt.scatter(pos_x_others_new, u.select_atoms('resname DSPC or resname CHL').positions[:,1])
#     plt.scatter(pos_x_DAPC_new, u.select_atoms('resname DAPC').positions[:,1])    
#     plt.show()
#     
# =============================================================================
    #plt.figure()
    #plt.scatter(u.select_atoms('resname DSPC or resname CHL').positions[:,0], u.select_atoms('resname DSPC or resname CHL').positions[:,1])
    #plt.scatter(u.select_atoms('resname DAPC').positions[:,0], u.select_atoms('resname DAPC').positions[:,1])    
    #plt.show()

# =============================================================================
# PI= math.pi
# radius_i = Lx/ (2. * PI)
# x_cm_i = np.mean(radius_i * np.cos(pos_DAPC_xy[:, 0] / Lx * 2. * PI))
# z_cm_i = np.mean(radius_i * np.sin(pos_DAPC_xy[:, 0] / Lx * 2. * PI))
# plt.scatter(x_cm_i, z_cm_i)
# print((np.arctan2(-z_cm_i, -x_cm_i)+PI)/ (2. * PI) * Lx)
# plt.axis('equal')
# plt.show()
# =============================================================================
