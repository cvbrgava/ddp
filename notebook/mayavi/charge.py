from numpy import *
from enthought.mayavi.mlab import *
from enthought.mayavi.api import Engine
    
a = -10.
b = 10.
ds = 1.

x, y, z = mgrid[a:b:ds, a:b:ds, a:b:ds];

q = [-1.0, -1.0, 1.5, 2.0] # set the charges
qpos = [[0.5, 2.5, 2.5], # set the positions. Note: All positions are offset by 0.5 to avoid calculating charges on the grid itself (causes 0-division!)
        [1.5, -2.5, -1.5],
        [3.5, -5.5, 3.5],
        [-2.5, -6.5, 0.5]]

def recalculate():
    global Ex
    global Ey
    global Ez
    Ex = 0 * x; # create our value arrays, filled with zeros
    Ey = 0 * y;
    Ez = 0 * z;

    for i in range(len(q)): # calculate the charge field from each electric charge
        r = (x - qpos[i][0])**2 + (y - qpos[i][1])**2 + (z - qpos[i][2])**2
        Ex = Ex + (q[i] * (x - qpos[i][0])) / (r)**1.5
        Ey = Ey + (q[i] * (y - qpos[i][1])) / (r)**1.5
        Ez = Ez + (q[i] * (z - qpos[i][2])) / (r)**1.5

recalculate()

# Now, let's prepare some output

fig = figure(fgcolor=(1,1,1), bgcolor=(0,0,0)) # set the background and foreground of our scene
camera = fig.scene.camera
camera.yaw(120)

streams = []
          
for s in range(len(q)):
    fl = flow(x,y,z,Ex,Ey,Ez,seed_scale=0.3, seed_resolution=6); # create another flow object
    fl.seed.widget.center = qpos[s]; # make seed sphere surround the charge
    fl.stream_tracer.initial_integration_step = 0.04 # make small integration steps for smooth lines and to avoid ugly lines
    fl.stream_tracer.maximum_propagation = 10.0 # avoid ugly lines by limitiing propagation
    fl.stream_tracer.integration_direction = 'both' # integrate in both directions for best results
    fl.seed.widget.enabled = False # remove the widget
    fl.actor.property.opacity = 0.5 # set the line opacity to make it semi-transparent
    streams.append(fl);
fig.scene.disable_render = True # avoid drawing on screen while rendering
fig.scene.off_screen_rendering = True
# animate
frames = 300
azimuth = 0
da = 360. / frames
for fi in range(1,frames):
    print fi
    qpos[0][1] = qpos[0][1] - 0.021
    qpos[1][2] = qpos[1][2] - 0.031
    qpos[3][0] = qpos[3][0] + 0.031
    recalculate()
    for s in range(len(q)):
        streams[s].mlab_source.u = Ex
        streams[s].mlab_source.v = Ey
        streams[s].mlab_source.w = Ez
        streams[s].seed.widget.enabled = True # remove the widget
        streams[s].seed.widget.center = qpos[s]; # make seed sphere surround the charge
        streams[s].seed.widget.enabled = False # remove the widget
    #azimuth = azimuth + da
    view(focalpoint=(0,0,0),azimuth=azimuth, elevation=90.)
    fnum = '%04d' % fi
    savefig("anim/output" + fnum + ".png",(1280,720))
    
fig.scene.disable_render = False # end avoid drawing on screen while rendering
fig.scene.off_screen_rendering = False

print "Done."
