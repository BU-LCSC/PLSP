#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:27:44 2021

@author: milechin
"""

import sys
import os
import json
import requests
from requests.auth import HTTPBasicAuth
import pathlib
import time
import errno
import pandas


# Setup Planet Data API base URL
base_url = "https://api.planet.com/data/v1"
orders_url = 'https://api.planet.com/compute/ops/orders/v2'
stats_url = "{}/stats".format(base_url)
quick_url = "{}/quick-search".format(base_url)

def setup_filter(coords, minyear, maxyear):
    
    geometry_filter = {"type": "GeometryFilter",
        "field_name": "geometry",
        "config": {
          "type": "Polygon",
          "coordinates": coords
        }
    }
    
    date_filter = {
        "type": "DateRangeFilter", # Type of filter -> Date Range
        "field_name": "acquired", # The field to filter on: "acquired" -> Date on which the "image was taken"
        "config": {
            "gte": "{}-01-01T00:00:00.000Z".format(minyear), # "gte" -> Greater than or equal to
            "lt":"{}-01-01T00:00:00Z".format(maxyear)
            }
        }
    

    ground_control =  {
                    "type": "StringInFilter",
                    "config": ["true"],
                    "field_name": "ground_control"
                  }
    quality_category = {
                    "type": "StringInFilter",
                    "config": ["standard"],
                    "field_name": "quality_category"
                  }
    cloud_cover =  {
                     "type": "RangeFilter",
                     "field_name": "cloud_cover",
                     "config": {"gte": 0, "lte": 0.5}
                  }
    asset = {
            "type": "AssetFilter",
            "config": [
                "ortho_analytic_4b_sr", "ortho_analytic_4b", "ortho_udm2"
            ]
        }


    permission = {
        "type":"PermissionFilter",
        "config":[
            "assets:download"
            ]
        }

    and_filter = {
        "type": "AndFilter",
        "config": [ 
            cloud_cover,
            quality_category,
            ground_control,
            date_filter,
            geometry_filter,
            permission,
            asset
            ]
        }
    
    return and_filter

def read_geometry(path):
    
    if(False == os.path.exists(path)):
        sys.exit("GeoJSON path doesn't exist: {}".format(path))
        
    with open(path, "r") as file1:
        geo = json.load(file1)
        
    return geo
    

def place_order(request, auth, order_name):
    print('Placing order!')
    headers = {'content-type': 'application/json'}
    
    response = requests.post(orders_url, data=json.dumps(request), auth=auth, headers=headers)
    print(response)
    print(response.reason)
    
    if(response.status_code != 202):
        print("Order failed for {}".format(order_name))
        return -1
    order_id = response.json()['id']
    print(order_id)
    order_url = orders_url + '/' + order_id
    return order_url


def download_results(results, output_dir, order_name, overwrite=False):
    results_urls = [r['location'] for r in results]
    results_names = [r['name'] for r in results]
    print('{} items to download'.format(len(results_urls)))
    
    failed_count = 0
    skipped_count = 0
    success_count = 0
    
    filename = "failed_downloads_{}.txt".format(order_name)
    
    with open(os.path.join(output_dir, filename), "w") as file1:
    
        for url, name in zip(results_urls, results_names):
            path = pathlib.Path(os.path.join(output_dir, 'data', name))

            if overwrite or not path.exists():
                print('downloading {} to {}'.format(name, path))
                #logging.info('downloading {} to {}'.format(name,path))
                r = requests.get(url, allow_redirects=True)
                if(r.status_code == 200):
                    path.parent.mkdir(parents=True, exist_ok=True)
                    open(path, 'wb').write(r.content)
                    success_count += 1
                else:
                    #logging.error('Status code {}, {} not downloaded.')
                    print('Status code {}, {} not downloaded.'.format(r.status_code, name))
                    failed_count += 1
                    file1.write("{}, {}, {}, {} \n".format(time.strftime("%Y%m%d-%H%M%S"), r.status_code, name, url))
            else:
                print('{} already exists, skipping {}'.format(path, name))
                skipped_count += 1
                #logging.info('{} already exists, skipping {}'.format(path, name))
                
    print("\n {} Success: {} Skipped: {} Failed: {}".format(order_name, success_count, skipped_count, failed_count))
    
    return failed_count, skipped_count, success_count

# Helper function to printformatted JSON using the json module
def p(data):
    print(json.dumps(data, indent=2))

def poll_for_success(order_url, auth):
    count = 0
    end_states = ['success', 'failed', 'partial']
    state = ''
    
    while(state not in end_states):
        count += 1
        time.sleep(10)
        r = requests.get(order_url, auth=auth)
        if(r.status_code != 200):
            print("\n Status code: {}".format(r.status_code))
            continue
        
        response = r.json()
        state = response['state']
        print(state)
        if state in end_states:
            return r
            break


def main(argv):
    start_time = time.time()
    
    site_num = pandas.to_numeric(argv[0])
    geometry_path = "/projectnb/modislc/users/seamorez/HLS_FCover/PLSP/geojson"
    min_year=argv[1]
    max_year=argv[2]
    #min_year = 2022
    #max_year = 2023
    output_dir = "/projectnb/modislc/users/seamorez/HLS_FCover/PLSP/rawImage"
    check_existing_orders="True"
    
    print("Input Params:")
    print("geometry_path={}".format(geometry_path))
    print("min_year={}".format(min_year))
    print("max_year={}".format(max_year))
    print("output_dir={}".format(output_dir))
    print("check_existing_orders={}".format(check_existing_orders))
    
    if(os.path.isdir(output_dir) == False):
        sys.exit("{} does not exist.".format(output_dir))
    
    #PLANET_API_KEY = os.getenv('PL_API_KEY')
    #PLANET_API_KEY = "3feab9a6cc1b4c1e8d65281025ad3382"
    PLANET_API_KEY = "PLAK8851ade7849744d294893d009ef4b87f"
    
    # Setup the session
    session = requests.Session()
    # Authenticate
    session.auth = (PLANET_API_KEY, "")
    
    # Make a GET request to the Planet Data API
    res = session.get(base_url)
    # Response status code
    if(res.status_code != 200):
        print("Cannot cannot to base server {} with status code {}".format(base_url, res.status_code))
        sys.exit("Cannot cannot to base server {} with status code {}".format(base_url, res.status_code))
    else:
        print("Base server is alive.")
        
    
    paths_list = os.listdir(geometry_path)
    paths_list.sort()
    print("\n Number of sites: {}".format(len(paths_list)))
    #print("Processing: {}".format(len(paths_list[10:])))
    #paths_list = paths_list[17:20]   
    #paths_list = [paths_list[i] for i in [37,38,39,40]]
    
    #for file in paths_list:
    file = paths_list[site_num]
    if file.endswith(".geojson"):
            
            print("\n Processing file: {}".format(file))
            
            file_path = os.path.join(geometry_path, file)
            output_site_dir = os.path.join(output_dir, os.path.splitext(file)[0]) 
            
            # if file == "Harvard_Forest_Hemlock_Site.geojson":
            #     check_existing_orders="True"
            # else:
            #     check_existing_orders="False"
            
            if not os.path.exists(output_site_dir):
                try:
                    print("\n Make directory")
                    os.makedirs(output_site_dir)
                except OSError as exc: # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            
            print("Output dir updated to: {}".format(output_site_dir))
            
            
            geo = read_geometry(file_path)
            
            for x in geo['features']:
                feature_name = x['properties']['f']
                feature_coords = x['geometry']['coordinates']
                
                filter = setup_filter(feature_coords, min_year, max_year)
                
                print("Filter config for search:")
                p(filter)
                
                print("\n\nRequesting stat count.")
                request = {
                    "interval" : "year",
                    "filter" : filter,
                    "item_types" : ["PSScene"]
                    }
        
                # Send the POST request to the API stats endpoint
                res = session.post(stats_url, json=request)
                
                if(res.status_code != 200):
                    sys.exit("Stats search failed with code {}".format(res.status_code))
        
                # Print response
                for bucket in res.json()['buckets']:
                    print("start_time: {} count: {}".format(bucket["start_time"], bucket["count"]))
                    
         
                print("do a real search.")
                request = {
                    "filter" : filter,
                    "item_types" : ["PSScene"]
                    }
        
                # Send the POST request to the API quick search endpoint
                res = session.post(quick_url, json=request)
                
                if(res.status_code != 200):
                    sys.exit("Quick search failed with code {}".format(res.status_code))
                    
                geojson = res.json()
                
                print("\n Creating file to save quick search results.")
                filename = "{}_quick_search_result_{}_{}.json".format(feature_name.replace(" ", "_"), min_year, max_year)
                with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                    json.dump(geojson, outfile)
                    
                
                print("\n Assemble ID list for download")
                features = geojson["features"]
                
                if(len(features) == 0):
                    sys.exit("0 IDs returned in quick search.")
                id_list = []
                while(len(geojson["features"]) > 0):
                    
                    for x in geojson["features"]:
                        id_list.append(x["id"])
                    
                    # Assign the "_links" -> "_next" property (link to next page of results) to a variable 
                    next_url = geojson["_links"]["_next"]
                    
                    time.sleep(5)
                    res = session.get(next_url)
                    
                    if(res.status_code != 200):
                        sys.exit("Next page retrieval failed with code {}".format(res.status_code))
                        
                    geojson = res.json()
                    with open(os.path.join(output_site_dir, filename), 'a') as outfile:
                        json.dump(geojson, outfile)
                
                print("\n Total IDs: {}".format(len(id_list)))
                
                print("\n Chunking ID list.")
                chunks = [id_list[x:x+400] for x in range(0, len(id_list), 400)]
                
                print("\n Number of chunks: {}".format(len(chunks)))
                
                if(len(chunks) >= 80):
                    sys.exit("{} Chunks which is greater than 80.  This will exceed order capacity".format(len(chunks)))
                
                #continue    
                #sys.exit()
            
                print("\n Checking connection with order server...")
                auth = HTTPBasicAuth(PLANET_API_KEY, '')
                response = requests.get(orders_url, auth=auth)
                
                if(response.status_code != 200):
                    sys.exit("Failed to connect to order server with code {}".format(response.status_code))
                else:
                    print("connected!")
                    
                orders_list = response.json()["orders"]
                    
                if(check_existing_orders == "True"):
                    
                    
                    while("next" in response.json()["_links"]):
                        
                        time.sleep(5)
                        next_url = response.json()["_links"]["next"]
                        response = session.get(next_url)
                        print(response.status_code)
                        if(response.status_code == 200):
                            
                            if("orders" in response.json()):
                            #print(len(response.json()["orders"]))
                                orders_list.extend(response.json()["orders"])
                        
                        print(len(orders_list))
                
                print("\n Time to place order...")
                
                orders_url_list = []
                bad_order_count = 0
                
                for count, chunk in enumerate(chunks):
                    print("\n Processing chunk: {}".format(count))
                    
                    order_name = "{}_chunk_{}_{}_{}".format(feature_name.replace(" ", "_"), count, min_year, max_year)
                    
                    
                    if(check_existing_orders == "True"):
                        
                        print("\n Checking existing orders")
                        orders_list = response.json()["orders"]
                        
                        for order in orders_list:
                            if(order["name"] == order_name):
                                print("\n Found existing Order by name.")
                                orders_url_list.append(order["_links"]["_self"])
                        
                    
                        feature_coords_buffer = feature_coords
                        feature_coords_buffer[0][0][0] = feature_coords_buffer[0][0][0] - 0.0015
                        feature_coords_buffer[0][0][1] = feature_coords_buffer[0][0][1] - 0.0015    
                        feature_coords_buffer[0][1][0] = feature_coords_buffer[0][1][0] + 0.0015
                        feature_coords_buffer[0][1][1] = feature_coords_buffer[0][1][1] - 0.0015
                        feature_coords_buffer[0][2][0] = feature_coords_buffer[0][2][0] + 0.0015
                        feature_coords_buffer[0][2][1] = feature_coords_buffer[0][2][1] + 0.0015
                        feature_coords_buffer[0][3][0] = feature_coords_buffer[0][3][0] - 0.0015
                        feature_coords_buffer[0][3][1] = feature_coords_buffer[0][3][1] + 0.0015
                        feature_coords_buffer[0][4][0] = feature_coords_buffer[0][4][0] - 0.0015
                        feature_coords_buffer[0][4][1] = feature_coords_buffer[0][4][1] - 0.0015
                        
                        request = {  
                           "name": order_name,
                           "order_type": "partial",
                           "products":[
                              {  
                                 "item_ids":chunk,
                                 "item_type":"PSScene",
                                  
                                 "product_bundle":"analytic_sr_udm2"
                              }
                           ],
                           "tools": [
                               {
                                 "clip": {     
                                     "aoi": {
                                         "type": "Polygon",
                                         "coordinates": feature_coords_buffer
                                         }
                                     }
                                 }
                           ]
                        }
                        
                        filename = "order_{}.json".format(order_name)
                        with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                            json.dump(request, outfile)
                        
                        time.sleep(5)
                        
                        order_result = place_order(request, auth, order_name)
                        
                        if(order_result == -1):
                            bad_order_count = bad_order_count +1
                            continue
                        else:
                            print("\n Placed order. Order url: {}".format(order_result))
                            orders_url_list.append(order_result)
                    
                    
                       
                        
                
                print("\n {} out of {} orders not placed successfully.".format(bad_order_count, len(chunks)))      
                
                if(len(orders_url_list) == 0):
                    sys.exit("No orders placed successfully.")
                
                    
                
                total_files_to_download = 0
                total_failed = 0
                total_skipped = 0
                total_success = 0
                order_list_failed = []
                order_list_failed_download = []
                
                print("\n Monitor orders...")
                for num, order in enumerate(orders_url_list):
                    try:
                        order_name = "{}_chunk_{}_{}_{}".format(feature_name.replace(" ", "_"), num, min_year, max_year)
                        
                        print("Monioring order num {}, url: {}".format(num, order))
                     
                        r = poll_for_success(order, auth)
                        
                        response = r.json()
                        state = response['state']
                        
                        if(state != "success"):
                            print("\n Order not success, but status {}".format(state))
                            order_list_failed.append(num)
                            continue
            
                        results = response['_links']['results']
                        
                        print("\n Total files to download: {} ".format(len(results)))
                        
                        total_files_to_download += len(results)
                        
                        filename = "order_result_{}.json".format(order_name)
                        with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                            json.dump(response, outfile)
            
                        failed_count, skipped_count, success_count = download_results(results, output_site_dir, order_name)
                        
                        if(failed_count != 0):
                            order_list_failed_download.append(num)
                        
                        total_failed += failed_count
                        total_skipped += skipped_count
                        total_success += success_count
                
                    except :
                        pass
                
                print("\nSUMMARY")
                print("Total Files to be downloaded: {}".format(total_files_to_download))
                print("Total Files Failed to Download: {}".format(total_failed))
                print("Total Files Skipped: {}".format(total_skipped))
                print("Total Files Success: {}".format(total_success))
                print("Order Chunks Failed: {}".format(order_list_failed))
                print("Download Chunks Failed: {}".format(order_list_failed_download))
                print("--- %s seconds ---" % (time.time() - start_time))
    
    else:
      print('No geojsons found!')
        
    
print(sys.argv)

if __name__ == "__main__":
    
   main(sys.argv[1:])

   
