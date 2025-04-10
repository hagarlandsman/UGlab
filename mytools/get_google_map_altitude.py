import requests
import json
from math import radians, sin, cos, acos, pi
from global_land_mask import globe

def jprint(obj):
    # create a formatted string of the Python JSON object
    text = json.dumps(obj, sort_keys=True, indent=4)
    print (text)
def get_alt(lat0, lon0, lat1,lon1,Nsteps):
    try:
       # response = requests.get("https://maps.googleapis.com/maps/api/elevation/json?locations="+str(lat)+"%2C"+str(long)+"&key=AIzaSyCtMSbCJeALFxBPlCNvB4Ruk_oVq0eBjKM")
        response = requests.get("https://maps.googleapis.com/maps/api/elevation/json?path="+str(lat0)+"%2C"+str(lon0)+"|"+str(lat1)+"%2C"+str(lon1)+"&samples="+str(Nsteps)+"&key=AIzaSyCtMSbCJeALFxBPlCNvB4Ruk_oVq0eBjKM")
        if (response.status_code == 200):
            j = response.json()
            return j
        else:
            return response.status_code
    except Exception as e:
        print ("problem ",lat0, lon0, lat1,lon1, e)
        return -1

def get_area(scan_radius, res, do_land =0, fout_name="", Nmax=500):
    if (fout_name==""):
        fout_name="R"+str(scan_radius)+"km_"+str(res)+".csv"
    fout = open(fout_name,"a")
    fout.write("latitude,longitude,altitude,x,y,land\n")
    fout.close()
    lat_lab = 32.60199361926758  #33.162712565    # power room
    lon_lab = 35.53016839929137  #35.525451159    #power room
    lab_alt = -280

    Rlat = 110.574;
    Rlon = 111.320 * cos(radians(lat_lab)) #0.837;
   # Clati = LatKochav * Rlati;
   # Clong = LonKochav * Rlong;


    lat_step = res/Rlat
    lon_step = res/Rlon

    lat_radius = scan_radius/Rlat
    lon_radius = scan_radius/Rlon

    lat_N = 2*lat_radius / lat_step
    lon_N = 2*lon_radius / lon_step

    print ("For range %f km, in res of %f km, lat_range is %f, step %f, N=%f "%(scan_radius,res,lat_radius,lat_step,lat_N))
    print ("For range %f km, in res of %f km, lon_range is %f, step %f, N=%f "%(scan_radius,res,lon_radius,lon_step,lon_N))


    lon_ranges=[]
    lon_this = lon_lab - lon_radius
    # Since up to Nmax points are allowed in each query, we will split the range to several queries.
    # the following loop puts into list lon_ranges the lower and upper bound of each segment.
    # while lon_this < lon_lab + lon_radius:
    while 1 :
        if (lon_lab + lon_radius < lon_this + Nmax*lon_step):
            lon_ranges.append((lon_this, lon_lab + lon_radius))
            break
        else:
            lon_ranges.append((lon_this, lon_this + Nmax*lon_step))
            lon_this = lon_this + (Nmax+1)*lon_step
        #lon_ranges.append((lon_this,min(lon_this + Nmax*lon_step, lon_this + lon_radius)))
        #lon_this = lon_this + lon_step*(Nmax+1)
    print (lon_ranges)
    lat_this = lat_lab - lat_radius
    l=[]
    c=0
    isl=-1
    while lat_this <= lat_lab+lat_radius:
        c=c+1
        l=[]
        for (lon0,lon1) in lon_ranges:
            ll = (get_alt(lat_this, lon0, lat_this, lon1, int((lon1-lon0)/lon_step)))
    #        print (ll)
            l.append(ll)
            print("%d / %d \t get_alt(%f,%f,%f,%f,%f)"%(c,lon_N,lat_this, lon0, lat_this, lon1, int((lon1-lon0)/lon_step)))
        with open(fout_name,"a") as fout:
                for eleres in l:
                    if (eleres['status']=='OK'):
                        for item in eleres['results']:
                            elevation = item['elevation']
                            llat = item['location']['lat']
                            lng = item['location']['lng']
                            if do_land == 1:
                                if globe.is_land(llat,lng):
                                    isl = 1
                                else:
                                    isl = 0

                            fout.write("%f,%f,%f,%f,%f,%d\n"%(llat,lng,elevation,(llat-lat_lab)*Rlat,(lng-lon_lab)*Rlon,isl))
                            #print ("%f %f %f" % ((llat-lat_lab)*Rlat,(lng-lon_lab)*Rlon,elevation))
        lat_this = lat_this + lat_step

if __name__ == '__main__':
#    get_area(scan_radius=10, res=0.1,fout_name="aa.csv", Nmax=10)
    get_area(scan_radius=100, res=1)
