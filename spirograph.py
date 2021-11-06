import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

plt.style.use("dark_background")

fig, ax = plt.subplots(figsize=(10, 5))
plt.axis("off")

ax.set_aspect("equal")
ax.set_xlim([-500, 500])
ax.set_ylim([-500, 500])

# Radius of the first circle, scales the whole image.
R = 40
# Ratio of the radius of the first circle to the second
k = 1
# Ratio of the radius of the first circle to the third
k2 = 2.01
# Distance to the drawing point
h = 40
# Ratio of the paths
p = 0.5

r = R / k
r2 = R / k2


# Length of the spiral
# High length can be difficult for the processor
len = 200
# If you see dots, decrease delta.
# Low delta can be difficult for the processor
delta = 0.001

tRange = np.arange(0, len, delta)


def hipo_draw(k, k2, h, p, t):
    """
    Draws the hypo(cycloid/trochoid) for two circles
    :param k: The ratio of the first (innermost) circle's radius to the second circle's radius
    :param k2: The ratio of the first circle's radius to the third (outermost) circle's radius
    :param h: The distance from the center of the third (outermost) circle to the drawing point
    :param p: The ratio of the arc length for rotation alpha of the innermost circle to the travelled arc length of the outermost circle
    :param t: The parametric range to draw the graph for
    :returns: The x and y coordinates for the graph drawn using the specified parameters
    """
    r = R / k
    x = (R - r) * np.cos(t) + h * np.cos(t - R * t / r)
    y = (R - r) * np.sin(t) + h * np.sin(t - R * t / r)
    return x, y


def hipo2_draw(k, k2, h, p, t):
    """
    Draws the hypo(cycloid/trochoid) for three circles
    :param k: The ratio of the first (innermost) circle's radius to the second circle's radius
    :param k2: The ratio of the first circle's radius to the third (outermost) circle's radius
    :param h: The distance from the center of the third (outermost) circle to the drawing point
    :param p: The ratio of the arc length for rotation alpha of the innermost circle to the travelled arc length of the outermost circle
    :param t: The parametric range to draw the graph for
    :returns: The x and y coordinates for the graph drawn using the specified parameters
    """
    r = R / k
    r2 = r / k2
    x = (R + r) * np.cos(t) - (r - r2) * np.cos(-t) + h * np.cos(t + r * t / r2)
    y = (R + r) * np.sin(t) - (r - r2) * np.sin(-t) + h * np.sin(t + r * t / r2)
    print(R, r, r2)
    return x, y


def epi_draw(k, k2, h, p, t):
    """
    Draws the epi(cycloid/trochoid) for two circles
    :param k: The ratio of the first (innermost) circle's radius to the second circle's radius
    :param k2: The ratio of the first circle's radius to the third (outermost) circle's radius
    :param h: The distance from the center of the third (outermost) circle to the drawing point
    :param p: The ratio of the arc length for rotation alpha of the innermost circle to the travelled arc length of the outermost circle
    :param t: The parametric range to draw the graph for
    :returns: The x and y coordinates for the graph drawn using the specified parameters
    """
    r = R / k
    x = (R + r) * np.cos(t) - h * np.cos(t + R * t / r)
    y = (R + r) * np.sin(t) - h * np.sin(t + R * t / r)
    return x, y


def epi2_draw(k, k2, h, p, t):
    """
    Draws the epi(cycloid/trochoid) for three circles
    :param k: The ratio of the first (innermost) circle's radius to the second circle's radius
    :param k2: The ratio of the first circle's radius to the third (outermost) circle's radius
    :param h: The distance from the center of the third (outermost) circle to the drawing point
    :param p: The ratio of the arc length for rotation alpha of the innermost circle to the travelled arc length of the outermost circle
    :param t: The parametric range to draw the graph for
    :returns: The x and y coordinates for the graph drawn using the specified parameters
    """
    r = R / k
    r2 = R / k2

    x = (
        (R + r) * np.cos(t)
        + (r + r2) * np.cos(t + (R * t) / r - R * t / (p * r))
        + h * np.cos(t + R * t / r - R * t / (p * r) - R * t / (p * r2))
    )
    y = (
        (R + r) * np.sin(t)
        + (r + r2) * np.sin(t + (R * t) / r - R * t / (p * r))
        + h * np.sin(t + R * t / r - R * t / (p * r) - R * t / (p * r2))
    )
    return x, y


def update(val):
    """
    Function called on the update of any of the sliders, re-draws the image using the selected function
    :param val: The values of the sliders
    """
    x, y = epi2_draw(
        k_slider.val,
        k2_slider.val,
        h_slider.val,
        p_slider.val,
        np.arange(0, tMax_slider.val, delta),
    )

    sp2.set_xdata(x)
    sp2.set_ydata(y)
    fig.canvas.draw_idle()


x, y = epi2_draw(k, k2, h, p, tRange)
(sp2,) = ax.plot(x, y, linewidth=0.5, color="red", antialiased=True)


# Sliders
h_ax = plt.axes([0.25, 0.09, 0.65, 0.01])
h_slider = Slider(h_ax, "h", valmin=0, valmax=200.0, valinit=h, valstep=1)
h_slider.on_changed(update)

p_ax = plt.axes([0.25, 0.07, 0.65, 0.01])
p_slider = Slider(p_ax, "p", valmin=-5, valmax=10.0, valinit=p, valstep=0.01)
p_slider.on_changed(update)

k_ax = plt.axes([0.25, 0.05, 0.65, 0.01])
k_slider = Slider(k_ax, "k", valmin=0, valmax=10.0, valinit=k, valstep=0.01)
k_slider.on_changed(update)

k2_ax = plt.axes([0.25, 0.03, 0.65, 0.01])
k2_slider = Slider(k2_ax, "k2", valmin=0.1, valmax=10, valinit=k2, valstep=0.01)
k2_slider.on_changed(update)

tMax_ax = plt.axes([0.25, 0.01, 0.65, 0.01])
tMax_slider = Slider(tMax_ax, "t", valmin=0, valmax=500, valinit=tRange[-1], valstep=1)
tMax_slider.on_changed(update)

plt.show()
