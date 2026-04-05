#!/bin/bash
# Deployment script for genepicoloc Colocalization Viewer
#
# Usage: ./deploy.sh [OPTIONS]
#
# Options:
#   --help, -h      Show this help message
#   --dev           Dev mode: bind-mount code, restart without rebuild
#   --local         Full local build (self-contained image)
#   --status        Show status of running container
#   --stop          Stop the running container
#   --logs          Show container logs

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Settings
SHINY_PORT=3838
CONTAINER_NAME=genepicoloc-shiny
IMAGE_NAME=genepicoloc-shiny
DATA_PATH="${SCRIPT_DIR}/../../../data/atlas"

# Resolve to absolute path
DATA_PATH="$(cd "$DATA_PATH" 2>/dev/null && pwd)" || {
  echo -e "\033[0;31mError: Data path not found: ${DATA_PATH}\033[0m"
  exit 1
}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

show_help() {
  echo ""
  echo -e "${BLUE}genepicoloc Colocalization Viewer - Deployment Script${NC}"
  echo ""
  echo "Usage: ./deploy.sh [OPTIONS]"
  echo ""
  echo "Options:"
  echo "  --help, -h      Show this help message"
  echo "  --dev           Dev mode: bind-mount R code for fast iteration"
  echo "  --local         Full local build (self-contained image)"
  echo "  --status        Show status of running container"
  echo "  --stop          Stop the running container"
  echo "  --logs [N]      Show container logs (optional: last N lines)"
  echo ""
  echo "Workflow:"
  echo "  1. ./deploy.sh --dev       # Iterate on R code (restart, no rebuild)"
  echo "  2. ./deploy.sh --local     # Test full build before deploying"
  echo ""
}

show_status() {
  echo -e "${BLUE}=== Container Status ===${NC}"
  podman ps -a --filter "name=${CONTAINER_NAME}"
}

stop_container() {
  echo -e "${BLUE}=== Stopping Container ===${NC}"
  podman rm -f ${CONTAINER_NAME} 2>/dev/null || true
  echo -e "${GREEN}Container stopped${NC}"
}

show_logs() {
  local lines="${1:-100}"
  echo -e "${BLUE}=== Container Logs (last ${lines} lines) ===${NC}"
  podman logs ${CONTAINER_NAME} 2>&1 | tail -n "$lines"
}

deploy_dev() {
  echo -e "${BLUE}=== Dev mode: ${IMAGE_NAME} ===${NC}"
  echo "  Port:      127.0.0.1:${SHINY_PORT}"
  echo "  Container: ${CONTAINER_NAME}"
  echo "  Data:      ${DATA_PATH}"
  echo "  Code:      bind-mounted (edit and restart)"
  echo ""

  podman rm -f ${CONTAINER_NAME} 2>/dev/null || true

  # Build image if it doesn't exist
  if ! podman image inspect ${IMAGE_NAME} &>/dev/null; then
    echo "First run: building image..."
    podman build -t ${IMAGE_NAME} "$SCRIPT_DIR"
  fi

  podman run -d -p ${SHINY_PORT}:3838 \
    -v "${DATA_PATH}:/app/data:ro" \
    -v "${SCRIPT_DIR}/app.R:/app/app.R:ro" \
    -v "${SCRIPT_DIR}/R:/app/R:ro" \
    -v "${SCRIPT_DIR}/www:/app/www:ro" \
    --name ${CONTAINER_NAME} ${IMAGE_NAME}

  sleep 2
  echo ""
  echo -e "${GREEN}=== Dev mode running! ===${NC}"
  echo -e "App: ${BLUE}http://localhost:${SHINY_PORT}/${NC}"
  echo ""
  echo -e "${YELLOW}After editing R code:${NC}"
  echo "  podman restart ${CONTAINER_NAME}"
  echo ""
  echo -e "${YELLOW}After changing Dockerfile/packages:${NC}"
  echo "  podman rm -f ${CONTAINER_NAME} && ./deploy.sh --dev"
}

deploy_local() {
  echo -e "${BLUE}=== Full build: ${IMAGE_NAME} ===${NC}"
  echo "  Port:      127.0.0.1:${SHINY_PORT}"
  echo "  Container: ${CONTAINER_NAME}"
  echo "  Data:      ${DATA_PATH}"
  echo ""

  podman rm -f ${CONTAINER_NAME} 2>/dev/null || true

  podman build -t ${IMAGE_NAME} "$SCRIPT_DIR"

  podman run -d -p ${SHINY_PORT}:3838 \
    -v "${DATA_PATH}:/app/data:ro" \
    --name ${CONTAINER_NAME} ${IMAGE_NAME}

  sleep 3
  echo ""
  echo "Checking logs..."
  podman logs ${CONTAINER_NAME} 2>&1 | tail -5

  echo ""
  echo -e "${GREEN}=== Deployment complete! ===${NC}"
  echo -e "App: ${BLUE}http://localhost:${SHINY_PORT}/${NC}"
  echo ""
  echo "Useful commands:"
  echo "  ./deploy.sh --logs     Show container logs"
  echo "  ./deploy.sh --status   Check container status"
  echo "  ./deploy.sh --stop     Stop the container"
}

case "${1:-}" in
  --help|-h|"")
    show_help
    exit 0
    ;;
  --dev)
    deploy_dev
    ;;
  --local)
    deploy_local
    ;;
  --status)
    show_status
    ;;
  --stop)
    stop_container
    ;;
  --logs)
    show_logs "${2:-100}"
    ;;
  *)
    echo -e "${RED}Unknown option: $1${NC}"
    echo "Use --help to see available options"
    exit 1
    ;;
esac
