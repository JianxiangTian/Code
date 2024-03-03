from PIL import Image
import numpy as np

# Read in an image
image_path = './configuration.png'
image = Image.open(image_path)

# Adjusting image size
image = image.resize((150, 150))

# Convert to grayscale image
image = image.convert("L")

# Convert to NumPy array
image_array = np.array(image)

# Binarize the image
threshold = 128
binary_image = (image_array <= threshold).astype(np.uint8)

# Record the position of each black pixel
black_pixel_positions = np.argwhere(binary_image == 1)

# Store the results
output_file =  './poision.txt'
with open(output_file, 'w') as file:
    for pos in black_pixel_positions:
        file.write(f"{pos[0]}\t{pos[1]}\n")


