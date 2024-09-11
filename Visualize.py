import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# # Load the CSV data
# data = np.genfromtxt('grid_output.csv', delimiter=',', skip_header=1)

# # Extract the coordinates and values
# x = data[:, 0]
# y = data[:, 1]
# z = data[:, 2]
# values = data[:, 3]

# # Filter out the points with value 0 (if needed)
# # filtered_data = data[values > 0]
# # x = filtered_data[:, 0]
# # y = filtered_data[:, 1]
# # z = filtered_data[:, 2]
# # values = filtered_data[:, 3]

# # Create a 3D scatter plot
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')

# # Map the values to colors (1 and 2 to different colors)
# colors = np.where(values == 1, 'red', 'blue')

# # Plot the data
# scatter = ax.scatter(x, y, z, c=colors, marker='o', s=5)


# # Add labels and title
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('3D Grid Visualization')

# grid_flat = np.loadtxt('grid_output.csv', delimiter=',')

# # Define the dimensions of the 3D grid
# sizeX, sizeY, sizeZ = 50, 50, 50  # Replace with your grid dimensions
# grid = grid_flat.reshape((sizeX, sizeY, sizeZ))

# # Step 2: Select a specific slice
# slice_index = 10  # Replace with the slice index you want to visualize
# slice_data = grid[:, :, slice_index]

# # Step 3: Visualize the selected slice
# plt.imshow(slice_data, cmap='viridis', interpolation='none')
# plt.colorbar(label='Grid Value')
# plt.title(f'Slice {slice_index} of the 3D Grid')

# # Show the plot
# plt.show()

# Step 1: Read the CSV file using pandas
file_path = 'grid_output.csv'  # Replace with the path to your actual CSV file

# Read the CSV file, skipping the header (if any) or dealing with non-numeric data
df = pd.read_csv(file_path, delimiter=',')

# Convert the DataFrame to a numpy array
# If your data has headers or non-numeric columns, you may need to drop or select the right columns
grid_flat = df.to_numpy()

# Define the dimensions of the 3D grid
sizeX, sizeY, sizeZ = 500, 500, 50  # Replace with your grid dimensions

# Reshape the flat grid data into a 3D grid
grid = grid_flat.reshape((sizeX, sizeY, sizeZ))

# Step 2: Select a specific slice
slice_index = 10  # Replace with the slice index you want to visualize
slice_data = grid[:, :, slice_index]

# Step 3: Create a scatter plot
# Find the coordinates of the points with value 1 and 2
x, y = np.indices(slice_data.shape)

# Flatten the arrays for plotting
x_flat = x.flatten()
y_flat = y.flatten()
slice_data_flat = slice_data.flatten()

# Filter points based on value
red_points = slice_data_flat == 1
blue_points = slice_data_flat == 2

# Plot red points (value 1)
plt.scatter(x_flat[red_points], y_flat[red_points], color='red', label='1s', s=10)

# Plot blue points (value 2)
plt.scatter(x_flat[blue_points], y_flat[blue_points], color='blue', label='2s', s=10)

# Add title and legend
plt.title(f'Slice {slice_index} of the 3D Grid')
plt.legend()
plt.gca().invert_yaxis()  # Optional: invert Y-axis to match image display
plt.show()
