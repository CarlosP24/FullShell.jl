# TagBot Instructions

## Issue Identified

The TagBot workflow in this repository is currently **manually disabled**. This prevents TagBot from automatically creating GitHub releases when new versions are registered in the Julia General Registry.

## What Was Fixed

1. **Simplified permissions**: Reduced the permissions in `TagBot.yml` to only what's necessary:
   - `actions: read`
   - `checks: read`
   - `contents: write` (required to create tags and releases)
   - `issues: read` (required for issue comment triggers)
   - `pull-requests: read`
   - `statuses: read`

   Removed unnecessary permissions like `deployments`, `discussions`, `packages`, `pages`, `repository-projects`, and `security-events`.

2. **Configuration verified**: The workflow already has the recommended triggers:
   - `issue_comment` trigger (responds when JuliaTagBot comments on issues)
   - `workflow_dispatch` trigger (allows manual execution)

## How to Enable TagBot

Since the workflow is manually disabled, it **cannot be enabled through code changes**. The repository owner must enable it manually:

1. Go to https://github.com/CarlosP24/FullShell.jl/actions/workflows/TagBot.yml
2. Click the "..." menu (three dots) or the "Enable workflow" button
3. Confirm to enable the workflow

## How TagBot Works

Once enabled, TagBot will:
1. Listen for comments from JuliaTagBot on issues
2. When a new version is registered in the Julia General Registry, JuliaTagBot will comment on the TagBot trigger issue
3. The workflow will automatically create a corresponding Git tag and GitHub release
4. Documentation will be built for the new release (if configured)

## Testing TagBot

After enabling, you can test TagBot by:
1. Triggering it manually using the "Run workflow" button on the Actions page
2. Setting the `lookback` input to check for recent registrations (default is 3 days)

## Additional Notes

- The `ssh: ${{ secrets.DOCUMENTER_KEY }}` parameter is used for pushing documentation. If you don't have this secret set up or don't need it, TagBot will still work for creating tags and releases.
- Make sure the package is registered in the Julia General Registry for TagBot to function.
