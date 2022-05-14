# pill-detector
A tool to automatically detect whether the pill blister is empty or not. 

<div align="center">
    <img src="https://github.com/IsabellaPoles/pill-detector/blob/main/detect_result.png" width="400px"</img> 
</div>

A POLIMI team is trying to automatically analyze pictures from pill blisters to understand how many have been used. One pic (‘blister.jpg’) is given to be analyzed. 
The aims are the following:

- detect the external contours of the blister and compute its aspect ratio;
- remove the different illumination effect;
- search in each side for the expiration date and locate the side that contains it;
- detect the pills still in the blister and the missing ones.
