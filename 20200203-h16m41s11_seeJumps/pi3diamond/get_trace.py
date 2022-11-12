path = "D:\data\Charge_state"
import datetime
import time
pi3d.confocal.TraceLength = 4000
pi3d.confocal.ReadoutInterval = 0.01
pi3d.confocal.CountInterval = 0.0001
directory = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
os.makedirs(os.path.join(path,directory))
for i in range(3):
    time.sleep(pi3d.confocal.TraceLength*pi3d.confocal.ReadoutInterval)
    filename = 'trace' + str(i) + '.txt'
    fp = os.path.join(path,directory, filename)
    np.savetxt(fp,pi3d.confocal.C)
pi3d.confocal.TraceLength = 4000
