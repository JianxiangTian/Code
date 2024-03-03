import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# Read in coordinate data
url = '.\\results.txt'
x = []
y = []
f = open(url)
lines = f.readlines()
for line in lines:
    line = line.split('\t')
    x.append(float(line[0]))
    y.append(float(line[1]))
f.close()

# Set the radius R of the disk and the size of the image.
R = 0.03
size = 5
fig, ax = plt.subplots(figsize=(size, size))
ax.set_xlim(0-R, 1+R)
ax.set_ylim(0-R, 1+R)
ax.set_aspect('equal', adjustable='box')  # 使图形显示正方形

# Iterate through coordinate points and draw solid disks.
for i in range(len(x)):
    circle = plt.Circle((x[i], y[i]), R, color='black', alpha=1.0)  # alpha 设置透明度
    ax.add_patch(circle)

# Hide the axes
ax.axis('off')

# Save the image as a square image
plt.savefig('./configuration.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()