#!/bin/bash

echo "Docker Container and Data Cleanup Script"
echo "======================================="
echo "This script will:"
echo "1. List all running containers"
echo "2. Stop all running containers"
echo "3. Remove stopped containers"
echo "4. Remove unused volumes"
echo "5. Remove unused networks"
echo "6. Remove unused images"

# Function to prompt for confirmation
confirm() {
    read -r -p "${1:-Are you sure you want to proceed? [y/N]} " response
    case "$response" in
        [yY][eE][sS]|[yY])
            true
            ;;
        *)
            echo "Operation cancelled."
            exit 0
            ;;
    esac
}

# List all running containers
echo -e "\nCurrently running containers:"
docker ps

# Confirm before proceeding
confirm "Do you want to stop all running containers? [y/N] "

# Stop all running containers
echo "Stopping all running containers..."
docker stop $(docker ps -q)

# Remove all stopped containers
echo "Removing stopped containers..."
docker rm $(docker ps -a -q)

# Remove unused volumes
echo "Removing unused volumes..."
docker volume prune -f

# Remove unused networks
echo "Removing unused networks..."
docker network prune -f

# List all images
echo -e "\nExisting Docker images:"
docker images

# Confirm before removing images
confirm "Do you want to remove all unused Docker images? [y/N] "

# Remove unused images
echo "Removing unused images..."
docker image prune -a -f

echo -e "\nCleanup completed!"
echo "To verify:"
echo "- Running containers: $(docker ps -q | wc -l)"
echo "- Total containers: $(docker ps -a -q | wc -l)"
echo "- Images remaining: $(docker images -q | wc -l)"
