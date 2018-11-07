from hexgrid import *
import matplotlib.pyplot as plt
'''
IC: i.e. t = 0
Q_a = bathymetry
Q_th = 0; Q_th(j = source area) = some number
Q_v = 0; Q_v(source) = some number
Q_cj = 0; Q_cj(source) = some number
Q_cbj = fraction of each sediment present in bed
Q_d = thickness of soft sediment, which can be eroded
Q_o = 0
'''
x = 2
y = 2
Ny = Nx = 10
Nj = 1

Q_th  = np.zeros((Ny, Nx))  # Turbidity current thickness
Q_v   = np.zeros((Ny, Nx))  # Turbidity current speed (scalar)
Q_cj  = np.zeros((Ny, Nx, Nj))  # jth current sediment volume concentration
Q_cbj = np.zeros((Ny, Nx, Nj))  # jth bed sediment volume fraction
Q_d   = np.ones((Ny, Nx)) * np.inf  # Thickness of soft sediment
Q_d[1:-1, 1:-1] = 1.5
Q_o   = np.zeros((Ny, Nx, 6))  # Density current outflow

# Source area
Q_th[y,x] = 1.5
Q_v[y,x] = 0.2
Q_cj[y,x,0] = 0.3
Q_cbj[1:-1,1:-1] = 0.4
Q_d[1:-1,1:-1] = 1
Q_d[y,x] = 5

ICstates = [Q_th, Q_v, Q_cj, Q_cbj, Q_d, Q_o]

grid=Hexgrid(Ny,Nx,ICstates=ICstates,reposeAngle=np.deg2rad(30),terrain=None)


# # # grid.Q_o[y,x] =
# # # print(grid.NEIGHBORS[:,:,0])
# # print("self.Q_th =\n", grid.Q_th )
# print("self.Q_v  =\n", grid.Q_v  )
# # print("self.Q_cj =\n", grid.Q_cj )
# # print("self.Q_cbj=\n", grid.Q_cbj)
# print("self.Q_d  =\n", grid.Q_d  )
# print("self.Q_a  =\n", grid.Q_a  )
# # print("self.Q_o  =\n", grid.Q_o  )
def plotCA():
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, aspect='equal')
    points = ax.scatter(grid.X[:, :, 0].flatten(), grid.X[:, :, 1].flatten(), marker='h', c=grid.Q_a.flatten())
    fig.colorbar(points, fraction=0.026)
    ax.set_title('Terrain(x,y)')
    plt.ion()
    # plt.show()


for i in range(50):

    grid.time_step()
    # plotCA()



print("done!")
# print("After \n")
# # print("self.Q_th =\n", grid.Q_th )
# print("self.Q_v  =\n", grid.Q_v  )
# # print("self.Q_cj =\n", grid.Q_cj )
# # print("self.Q_cbj=\n", grid.Q_cbj)
# print("self.Q_d  =\n", grid.Q_d  )
# print("self.Q_a  =\n", grid.Q_a  )
# # print("self.Q_o  =\n", grid.Q_o  )
# print("dt=",grid.dt)
