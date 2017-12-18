import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from kelp import decompose

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
        
fig = plt.figure(figsize = (12,12))
ax = fig.gca(projection = '3d')
ax.set_aspect('equal')
fig.tight_layout()
MAX = 1
MIN = -.1

for point in np.diag(MAX * np.array([1,1,1])):
    ax.plot([point[0]], [point[1]], [point[2]], 'w')
    
for point in np.diag(MIN * np.array([1,1,1])):
    ax.plot([point[0]], [point[1]], [point[2]], 'w')


x, y = np.meshgrid(np.linspace(-.1,1,100),np.linspace(-.1,1,100))

z = (x+y)/2


a = np.array([6,2,4])/np.sqrt(54)
b = np.array([2,6,4])/np.sqrt(54)
c = np.array([2,3,7])/np.sqrt(62)

coeffs, r, chat = decompose(c,[a,b],full=True)

long_a_points = tuple(map(list,zip(np.array([0,0,0]),a*1.3)))
a_points = tuple(map(list,zip(np.array([0,0,0]),a)))
b_points = tuple(map(list,zip(np.array([0,0,0]),b)))
c_points = tuple(map(list,zip(np.array([0,0,0]),c)))
chat_points = tuple(map(list,zip(np.array([0,0,0]),chat)))
r_points = tuple(map(list,zip(c,chat)))

long_A = Arrow3D(*long_a_points, mutation_scale=10, lw=2, arrowstyle="-|>", color="g")
A = Arrow3D(*a_points, mutation_scale=10, lw=2, arrowstyle="-|>", color="g")
B = Arrow3D(*b_points, mutation_scale=10, lw=2, arrowstyle="-|>", color="b")
C = Arrow3D(*c_points, mutation_scale=10, lw=2, arrowstyle="-|>", color='r')
Chat = Arrow3D(*chat_points, mutation_scale=10, lw=2, ls='--', arrowstyle="-|>", color='r')
R = Arrow3D(*r_points, mutation_scale=10, lw=2, ls='--',arrowstyle="-|>", color='r')

ax.view_init(elev = 5, azim = 190)
ax.add_artist(long_A)
ax.add_artist(B)
# ax.set_xticks([])                               
# ax.set_yticks([])                               
# ax.set_zticks([])
plt.savefig('../plots/plot1.png')
long_A.remove()
ax.add_artist(A)
plt.savefig('../plots/plot2.png')
ax.plot_surface(x,y,z,shade=False,color='w',edgecolors='.75')
plt.savefig('../plots/plot3.png')
ax.add_artist(C)
plt.savefig('../plots/plot4.png')
ax.add_artist(Chat)
plt.savefig('../plots/plot5.png')
ax.add_artist(R)
plt.savefig('../plots/plot6.png')
R.remove()
C.remove()
plt.savefig('../plots/plot7.png')

# for ii in range(0,60,5):
    # ax.view_init(elev = ii, azim = 181)
    # plt.savefig('../plots/elev = ' + str(ii) + '.png')
# ax.set_xlim(0,10)
# ax.set_ylim(0,10)
# ax.set_zlim(0,10)
# plt.axis('equal')
plt.show()

