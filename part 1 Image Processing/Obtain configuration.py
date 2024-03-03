import cv2

# Read in an image
image = cv2.imread('./Image.jpg', cv2.IMREAD_GRAYSCALE)

# Obtain the width and height of the image
image_width = image.shape[1]
image_height = image.shape[0]
print("The width and height of the image are respectively","\t",image_height,"\t",image_width)

# Binarize the image
_, thresholded = cv2.threshold(image, 180, 255, cv2.THRESH_BINARY)

# Find each contour
contours, _ = cv2.findContours(thresholded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Store the results
with open('./results.txt', 'w') as file:
    for i, contour in enumerate(contours):
        # Calculate the geometric centroid of the contour
        M = cv2.moments(contour)
        if M["m00"] != 0:
            cx = int(M["m10"] / M["m00"])
            cy = int(M["m01"] / M["m00"])
        else:
            continue

        # Calculate the area and filter out factors with excessively small areas.
        area = cv2.contourArea(contour)
        if area > 10:
            # Normalize the coordinates
            cx = cx/image_width
            cy = cy/image_height
            file.write(f"{cx}\t{cy}\n")
