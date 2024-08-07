def build_correlation_matrix(M):
        '''
        Covert covariance matrix to correlation matrix

        Parameters:
        ----
        M (N pixels, k) : K = M M' , K is a covariance matrix, M is its matrix decomposition   

        Return:
        ----
        C : correlations matrix with its diag C = I
        '''
        # build covariance matrix
        K = np.matmul( M, M.T )

        # query diag elements
        d = np.sqrt(np.diag( K ))[:, np.newaxis]
        
        M_div_d = M / d

        C = np.matmul( M_div_d, M_div_d.T)

        return C

from scipy.io import loadmat                                                                            
import matplotlib.pyplot as plt
import numpy as np
M = loadmat("M.mat")                                                                                   
# M = loadmat("MM-70%-1350-1570.mat")                                                                                   
C = build_correlation_matrix(M["M"])    
# min_lambda = 1216
# max_lambda = 1600
# plt.imshow(C, origin='low')
# plt.show()
# # locs, labels= plt.xticks()
# # locs
# # new_labels = np.linspace(min_lambda, max_lambda,len(locs)-1)
# # new_labels = np.int16(new_labels)
# # print(new_labels)
# # x_ticks=[]
# # plt.xticks(locs, new_labels)
# # plt.savefig('MM-0.5-1350-1600.png')
# plt.savefig('MM-1216-1600.png')



# plotting covariance matrix
min_lambda = 1310;          
max_lambda = 1555
scale = np.shape(C)[0] / (max_lambda - min_lambda)
c2=1335
si4 = 1402
n4 = 1486
c4 = 1548
fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(C, origin="lower")
ax.set_xticks(
        [
            (c2 - min_lambda) * scale,
            (si4 - min_lambda) * scale,
            (n4 - min_lambda) * scale,
            (c4 - min_lambda) * scale,
        ],
    )
ax.set_xticklabels(
        [
            r"CII",
            r"SiIV",
            r"NIV",
            r"CIV",
        ],
        rotation=45,
    )
ax.set_yticks(
        [
            (c2 - min_lambda) * scale,
            (si4 - min_lambda) * scale,
            (n4 - min_lambda) * scale,
            (c4 - min_lambda) * scale,
        ]
    )
ax.set_yticklabels(
        [
             r"CII",
             r"SiIV",
            r"NIV",
            r"CIV",
        ]
    )
ax.tick_params(labelsize=40)
# fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
# plt.title('RATING<2')
plt.tight_layout()
# plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.savefig("covariance_matrix-M2")
plt.show()
