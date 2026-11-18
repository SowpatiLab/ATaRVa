from pathlib import Path

__version__ = "0.9.5"

def update_version(version):
    """
    Updating the version number in pyproject.toml and README.md files.

    :param version: The new version number.
    """
    lines = []

    path = str(Path(__file__).parent.resolve()) + '/../'
    with open(f"{path}pyproject.toml", "r") as f:
        for line in f:
            if line.startswith("version = "):
                line = f'version = "{version}"\n'
            lines.append(line)

    with open(f"{path}pyproject.toml", "w") as f:
        f.writelines(lines)

    lines = []
    with open(f"{path}lib/README.md", "r") as f:
        for line in f:
            if "{{VERSION}}" in line:
                line = line.replace("{{VERSION}}", version)
            lines.append(line)

    with open(f"{path}README.md", "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    print(f"Updating version to {__version__}")
    update_version(__version__)
