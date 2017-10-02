"""
========================================
An animated image using a list of images
========================================

This examples demonstrates how to animate an image from a list of images (or
Artists).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.close("all")

fig = plt.figure("Animation")
ax = fig.add_subplot(111)
ax.set_title('titre 1')


def f(x, y):
    return np.sin(x) + np.cos(y)

x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []
for i in range(60):
    x += np.pi / 15.
    y += np.pi / 20.
    img   = ax.imshow(f(x, y), animated=True)
    title = ax.annotate(i,(5,5)) # add text
    #title = ax.set_title("i = %2i" %(i))
    ims.append([img,title])

ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                repeat_delay=1000)

# ani.save('dynamic_images.mp4')

plt.show()

'''
ims = []
fig = plt.figure("Animation")
ax = fig.add_subplot(111)

for imgNum in range(10):
    img = np.random.rand(10,10) #random image for an example

    frame =  ax.imshow(img)
    t = ax.annotate(imgNum,(1,1)) # add text

    ims.append([frame,t]) # add both the image and the text to the list of artists

anim = animation.ArtistAnimation(fig, ims, interval=350, blit=True, repeat_delay=350)

plt.show()
'''