#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--publish] [--dry-run] [--no-changelog]

Options:
  --publish       Publish to crates.io after bumping and committing
  --dry-run       Show what would be done without modifying files, creating
                  commits or tags, or publishing
  --no-changelog  Skip regenerating CHANGELOG.md (requires git-cliff otherwise)
  -h, --help      Show this help message

CHANGELOG.md is regenerated from conventional commits by git-cliff (see
cliff.toml) and included in the release commit. Preview the section the next
release would add, without writing anything:

  git-cliff --tag v<version> --unreleased
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false
CHANGELOG_ENABLED=true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish)      PUBLISH=true ;;
        --dry-run)      DRY_RUN=true ;;
        --no-changelog) CHANGELOG_ENABLED=false ;;
        -h|--help)      usage; exit 0 ;;
        -*)             die "unknown option: $1" ;;
        *)
            [[ -z "$VERSION" ]] || die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || { usage; exit 1; }

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CRATE="grangers"
CARGO_TOML="Cargo.toml"
CHANGELOG="CHANGELOG.md"
CLIFF_CONFIG="cliff.toml"
TAG="v${VERSION}"

[[ -f "$CARGO_TOML" ]] || die "not found: $CARGO_TOML"

if [[ "$CHANGELOG_ENABLED" == true ]]; then
    [[ -f "$CLIFF_CONFIG" ]] || die "not found: $CLIFF_CONFIG (pass --no-changelog to skip changelog generation)"
    command -v git-cliff >/dev/null 2>&1 || die "git-cliff is not installed; install it (cargo binstall git-cliff) or pass --no-changelog"
fi

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$CARGO_TOML" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine current version from $CARGO_TOML"

[[ "$CURRENT_VERSION" != "$VERSION" ]] || die "crate version is already set to $VERSION"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi

if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi

echo "Crate            : $CRATE"
echo "Current version  : $CURRENT_VERSION"
echo "New version      : $VERSION"
echo "Tag              : $TAG"
echo "Publish crate    : $([[ "$PUBLISH" == true ]] && echo yes || echo no)"
echo "Regen $CHANGELOG : $([[ "$CHANGELOG_ENABLED" == true ]] && echo yes || echo no)"
echo "Dry-run          : $([[ "$DRY_RUN" == true ]] && echo yes || echo no)"
echo

echo "Updating $CARGO_TOML"
echo "  version: $CURRENT_VERSION -> $VERSION"

if [[ "$DRY_RUN" == false ]]; then
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$CARGO_TOML"
    rm -f "${CARGO_TOML}.bak"

    UPDATED_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$CARGO_TOML" | head -1)"
    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "${CARGO_TOML} version update failed"
else
    echo "Dry-run: would rewrite $CARGO_TOML"
fi

# Also refreshes Cargo.lock, which is not committed for this crate but must
# agree with the manifest for the checks below to mean anything.
run cargo check -q

# The parse-path fixtures are the only guard against a silent behaviour change
# in the GTF/GFF reader, so a release must not skip them.
run cargo test -q

if [[ "$CHANGELOG_ENABLED" == true ]]; then
    # Regenerate the whole file rather than prepending: the result is
    # idempotent and stays in one format throughout. `--tag` labels the
    # not-yet-tagged commits with the version about to be cut; without it they
    # would land under "Unreleased".
    echo "Regenerating $CHANGELOG for $TAG"
    run git-cliff --tag "$TAG" -o "$CHANGELOG"
fi

run git add "$CARGO_TOML"
if [[ "$CHANGELOG_ENABLED" == true ]]; then
    run git add "$CHANGELOG"
fi
run git commit -m "chore(release): bump ${CRATE} to v${VERSION}"

if [[ "$PUBLISH" == true ]]; then
    run cargo publish
fi

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$DRY_RUN" == true ]]; then
    echo
    echo "Dry-run complete"
else
    echo
    echo "Release complete for ${CRATE} v${VERSION}"
fi
