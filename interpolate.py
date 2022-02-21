import numpy as np
import plotly.graph_objects as go

def doInterpolation(data):
    for i in range(data.shape[0] - 1):
        data = interpolation(data[i,:], data[i+1,:], data)
    return data[data[:,3].argsort()]

def interpolation(p0, p1, data):
    x0, y0, z0, t0 = p0
    x1, y1, z1, t1 = p1
    unitVector = np.array([x1-x0, y1-y0, z1-z0])
    magnitudeUnit = np.linalg.norm(unitVector)
#     x = [x0]
#     y = [y0]
#     z = [z0]
    howMany = int(t1 - t0)
#     print(howMany)
    if (howMany > 1):
        for i in range(1, howMany):
            tempx = x0 + (magnitudeUnit*i/howMany)*(x1-x0)/(magnitudeUnit)
            tempy = y0 + (magnitudeUnit*i/howMany)*(y1-y0)/(magnitudeUnit)
            tempz = z0 + (magnitudeUnit*i/howMany)*(z1-z0)/(magnitudeUnit)
            tempt = t0 + i
            data = np.append(data, np.array([[tempx, tempy, tempz, tempt]]), axis=0)
#             print("Hello")
#     x.append(x1)
#     y.append(y1)
#     z.append(z1)
    return data

def plotCoord(x, y, z):
    fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z,
                                   mode='markers')])
    fig.show()

def interpole(p0, p1, i):
    x0, y0, z0, t0 = p0
    x1, y1, z1, t1 = p1
    unitVector = np.array([x1-x0, y1-y0, z1-z0])
    magnitudeUnit = np.linalg.norm(unitVector)

    howMany = int(t1 - t0)
    tempx = x0 + (magnitudeUnit*(i - t0)/howMany)*(x1-x0)/(magnitudeUnit)
    tempy = y0 + (magnitudeUnit*(i - t0)/howMany)*(y1-y0)/(magnitudeUnit)
    tempz = z0 + (magnitudeUnit*(i - t0)/howMany)*(z1-z0)/(magnitudeUnit)
    tempt = i

    return np.array([[tempx, tempy, tempz, tempt]])

def doInter(data, timesteps):
    j = 0
    for i, t in enumerate(timesteps):
        if (t != data[j][3]):
                temp = interpole(data[j-1], data[j], t)
                data = np.append(data, temp, axis = 0)
        if(t == data[j][3]):
            j = j+1
    return data