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

# Load per-machine config from .env (gitignored). See .env.example.
if [ -f "${SCRIPT_DIR}/.env" ]; then
  set -a
  # shellcheck disable=SC1091
  . "${SCRIPT_DIR}/.env"
  set +a
fi

# Settings
SHINY_PORT="${SHINY_PORT:-3838}"
CONTAINER_NAME="${CONTAINER_NAME:-genepicoloc-shiny}"
IMAGE_NAME="${IMAGE_NAME:-genepicoloc-shiny}"

# Atlas data path: GENEPICOLOC_DATA_PATH env var > repo-relative default.
# Set GENEPICOLOC_DATA_PATH in .env for local dev or in the systemd unit /
# docker-compose file for remote deploys.
DATA_PATH_RAW="${GENEPICOLOC_DATA_PATH:-${SCRIPT_DIR}/../../../data/atlas}"

# Resolve atlas path (only required for --dev / --local)
require_data_path() {
  DATA_PATH="$(cd "$DATA_PATH_RAW" 2>/dev/null && pwd)" || {
    echo -e "\033[0;31mError: atlas data path not found: ${DATA_PATH_RAW}\033[0m"
    echo "Set GENEPICOLOC_DATA_PATH in ${SCRIPT_DIR}/.env (see .env.example)"
    exit 1
  }
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
  echo "  --remote        Build source bundle for VPS (synced via Nextcloud)"
  echo "  --status        Show status of running container"
  echo "  --stop          Stop the running container"
  echo "  --logs [N]      Show container logs (optional: last N lines)"
  echo ""
  echo "Workflow:"
  echo "  1. ./deploy.sh --dev       # Iterate on R code (restart, no rebuild)"
  echo "  2. ./deploy.sh --local     # Test full build before deploying"
  echo "  3. ./deploy.sh --remote    # Bundle source -> Nextcloud -> VPS"
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
  require_data_path
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
    -e ANTHROPIC_API_KEY \
    -v "${DATA_PATH}:/app/data:ro" \
    -v "${SCRIPT_DIR}/app.R:/app/app.R:ro" \
    -v "${SCRIPT_DIR}/R:/app/R:ro" \
    -v "${SCRIPT_DIR}/www:/app/www:ro" \
    -v "${SCRIPT_DIR}/extdata:/app/extdata:ro" \
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
  require_data_path
  echo -e "${BLUE}=== Full build: ${IMAGE_NAME} ===${NC}"
  echo "  Port:      127.0.0.1:${SHINY_PORT}"
  echo "  Container: ${CONTAINER_NAME}"
  echo "  Data:      ${DATA_PATH}"
  echo ""

  podman rm -f ${CONTAINER_NAME} 2>/dev/null || true

  podman build -t ${IMAGE_NAME} "$SCRIPT_DIR"

  podman run -d -p ${SHINY_PORT}:3838 \
    -e ANTHROPIC_API_KEY \
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

# Remote deployment: build a source-only tarball that the VPS will rebuild.
# Mirrors apps/nako/deploy.sh + deploy-remote.sh pattern.
create_bundle() {
  BUNDLE_NAME="genepicoloc-shiny.tar.gz"
  BUNDLE_PATH="${SCRIPT_DIR}/${BUNDLE_NAME}"

  echo -e "${BLUE}=== Creating source bundle ===${NC}"
  cd "$SCRIPT_DIR"
  tar -czf "$BUNDLE_PATH" \
    --exclude='*.tar.gz' \
    --exclude='.git' \
    --exclude='.Rhistory' \
    --exclude='.env' \
    Dockerfile \
    app.R \
    R \
    www \
    extdata \
    .env.example

  BUNDLE_SIZE=$(du -h "$BUNDLE_PATH" | cut -f1)

  # Drop into the Nextcloud-synced epi-vps tree so the VPS sees it.
  VPS_APPS_DIR="${HOME}/Work/bioinfo/nextcloud/Downloads/epi-vps/apps/genepicoloc"
  if [[ -d "$VPS_APPS_DIR" ]]; then
    mv "$BUNDLE_PATH" "$VPS_APPS_DIR/"
    echo -e "${GREEN}Bundle copied to: ${VPS_APPS_DIR}/${BUNDLE_NAME} (${BUNDLE_SIZE})${NC}"
  else
    echo -e "${GREEN}Bundle created: ${BUNDLE_PATH} (${BUNDLE_SIZE})${NC}"
    echo -e "${YELLOW}VPS sync folder not found: ${VPS_APPS_DIR}${NC}"
  fi

  echo ""
  echo "Next steps:"
  echo "  1. On laptop: bash ~/epi-vps/apps/genepicoloc/deploy.sh         (scp bundle to epi)"
  echo "  2. On epi:    proxy_on && bash ~/epi-vps/apps/genepicoloc/deploy-remote.sh"
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
  --remote|--bundle)
    create_bundle
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
