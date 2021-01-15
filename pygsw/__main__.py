import pygsw

if __name__ == '__main__':
    locations = ([0., 0., -5000, 10, 30],)
    for long, lat, z, t, sp in locations:
        pt = pygsw.calculate_pt(long, lat, z, t, sp)
        print(pt)

