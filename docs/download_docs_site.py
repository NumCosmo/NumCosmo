#!/usr/bin/env python3
"""
Download pre-built documentation site artifact from GitHub Actions.

This script is designed to run during Read the Docs builds. It downloads
the pre-built documentation site from GitHub Actions artifacts,
which significantly speeds up the documentation deployment by reusing
the site built during the GitHub Actions workflow.

Required environment variables:
    READTHEDOCS_GIT_COMMIT_HASH: The Git commit SHA for the current build
    GITHUB_TOKEN: GitHub token for API authentication (optional but recommended)
"""

import os
import shutil
import sys
import zipfile

import requests

# Configuration constants
REPO = "NumCosmo/NumCosmo"
SITE_ARTIFACT_PREFIX = "numcosmo-site"
API_TIMEOUT_SECONDS = 120
DOWNLOAD_TIMEOUT_SECONDS = 3600
MAX_PAGES = 50
SITE_ZIP_NAME = "site.zip"
SITE_TARGET_DIR = "numcosmo-site"


def http_get_json(url: str, headers: dict[str, str], timeout: int) -> dict:
    """
    Perform HTTP GET request and parse JSON response.

    Args:
        url: URL to fetch
        headers: HTTP headers for the request
        timeout: Request timeout in seconds

    Returns:
        Parsed JSON response as dictionary

    Raises:
        SystemExit: If request fails or response is not valid JSON
    """
    try:
        response = requests.get(url, headers=headers, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {url}: {e}")
        sys.exit(1)


def http_download_file(
    url: str,
    filename: str,
    headers: dict[str, str],
    timeout: int,
    chunk_size: int = 8192,
) -> None:
    """
    Download file from URL using streaming to handle large files.

    Args:
        url: URL to download from
        filename: Local filename to save to
        headers: HTTP headers for the request
        timeout: Request timeout in seconds
        chunk_size: Size of chunks for streaming download

    Raises:
        SystemExit: If download or file writing fails
    """
    try:
        with requests.get(
            url, headers=headers, timeout=timeout, stream=True
        ) as response:
            response.raise_for_status()
            with open(filename, "wb") as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        sys.exit(1)
    except IOError as e:
        print(f"Error writing file {filename}: {e}")
        sys.exit(1)


def get_configuration() -> tuple[str, str | None, dict[str, str]]:
    """
    Get configuration from environment variables and prepare API headers.

    Returns:
        Tuple of (artifact_name, token, headers)

    Raises:
        SystemExit: If required environment variables are missing
    """
    sha = os.environ.get("READTHEDOCS_GIT_COMMIT_HASH")
    token = os.environ.get("GITHUB_TOKEN")

    if not sha:
        print("Error: READTHEDOCS_GIT_COMMIT_HASH environment variable is required")
        sys.exit(1)

    artifact_name = f"{SITE_ARTIFACT_PREFIX}-{sha}"

    print(f"Looking for artifact: {artifact_name}")
    if token:
        print("Using authenticated GitHub API requests")
    else:
        print(
            "Warning: No GITHUB_TOKEN set, using "
            "unauthenticated requests (rate-limited)"
        )

    headers = {"Accept": "application/vnd.github+json"}
    if token:
        # Use 'token' format for GitHub PATs (not 'Bearer')
        # This format works for artifact downloads and API calls
        headers["Authorization"] = f"token {token}"

    return artifact_name, token, headers


def find_artifact(repo: str, artifact_name: str, headers: dict[str, str]) -> str | None:
    """
    Find the artifact by name using GitHub API with pagination.

    Args:
        repo: Repository in "owner/name" format
        artifact_name: Name of the artifact to find
        headers: HTTP headers for API requests

    Returns:
        Download URL of the artifact, or None if not found

    Raises:
        SystemExit: If API request fails
    """
    url = f"https://api.github.com/repos/{repo}/actions/artifacts"
    page = 1

    while page <= MAX_PAGES:
        paged_url = f"{url}?per_page=100&page={page}"
        print(f"Fetching artifacts page {page}")

        data = http_get_json(paged_url, headers, API_TIMEOUT_SECONDS)
        artifacts = data.get("artifacts", [])

        if not artifacts:
            break

        for artifact in artifacts:
            if artifact["name"] == artifact_name and not artifact["expired"]:
                return artifact["archive_download_url"]

        page += 1
    else:
        print(f"Reached maximum page limit ({MAX_PAGES}) without finding artifact")

    return None


def download_artifact(
    download_url: str, filename: str, headers: dict[str, str]
) -> None:
    """
    Download artifact from GitHub.

    Args:
        download_url: URL to download the artifact from
        filename: Local filename to save the artifact to
        headers: HTTP headers for the request

    Raises:
        SystemExit: If download or file writing fails
    """
    print(f"Downloading artifact from: {download_url}")
    http_download_file(download_url, filename, headers, DOWNLOAD_TIMEOUT_SECONDS)


def is_safe_zip_path(target_dir: str, member_path: str) -> bool:
    """
    Validate that a zip member path is safe to extract.

    Protects against Zip Slip attacks by rejecting:
    - Absolute paths
    - Paths containing .. segments
    - Paths that would resolve outside the target directory

    Args:
        target_dir: The intended extraction directory
        member_path: The path from the zip member

    Returns:
        True if the path is safe to extract, False otherwise
    """
    # Reject absolute paths
    if os.path.isabs(member_path):
        return False

    # Reject paths with .. segments
    if ".." in os.path.normpath(member_path).split(os.sep):
        return False

    # Verify the resolved path stays within target directory
    target_dir_abs = os.path.abspath(target_dir)
    member_abs = os.path.abspath(os.path.join(target_dir, member_path))

    # Check if the member path is within target directory
    return (
        member_abs.startswith(target_dir_abs + os.sep) or member_abs == target_dir_abs
    )


def extract_artifact(zip_filename: str, target_dir: str) -> None:
    """
    Safely extract artifact zip file to target directory and clean up.

    Validates all member paths to prevent Zip Slip attacks.

    Args:
        zip_filename: Path to the zip file
        target_dir: Directory to extract contents to

    Raises:
        SystemExit: If extraction fails or malicious paths are detected
    """
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)

    os.makedirs(target_dir, exist_ok=True)

    print(f"Extracting artifact to {target_dir}")
    try:
        with zipfile.ZipFile(zip_filename) as z:
            # Validate all member paths before extraction
            for member in z.namelist():
                if not is_safe_zip_path(target_dir, member):
                    print(f"Error: Unsafe path in zip file: {member}")
                    print("Possible Zip Slip attack detected")
                    sys.exit(1)

            # Safe to extract after validation
            z.extractall(target_dir)

    except (zipfile.BadZipFile, OSError) as e:
        print(f"Error extracting artifact: {e}")
        sys.exit(1)
    finally:
        # Clean up the downloaded zip file
        if os.path.exists(zip_filename):
            os.remove(zip_filename)
            print(f"Cleaned up {zip_filename}")


def main() -> None:
    """Main execution function."""
    # Get configuration
    artifact_name, _token, headers = get_configuration()

    # Find artifact
    download_url = find_artifact(REPO, artifact_name, headers)
    if download_url is None:
        print(f"Artifact not found: {artifact_name}")
        sys.exit(1)

    # Download artifact
    download_artifact(download_url, SITE_ZIP_NAME, headers)

    # Extract artifact
    extract_artifact(SITE_ZIP_NAME, SITE_TARGET_DIR)

    print(f"Successfully extracted documentation site into {SITE_TARGET_DIR}")


if __name__ == "__main__":
    main()
