# thebustalab.github.io

Install the pre-push secret scan once per clone (blocks a push whose outgoing diff
carries a credential or a private path, and names each hit as `path:line`):

```bash
ln -sf ../../githooks/pre-push .git/hooks/pre-push
```
