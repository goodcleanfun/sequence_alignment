
test:
	clib install --dev
	@$(CC) test.c src/local_alignment.c deps/utf8/utf8.c deps/utf8proc/utf8proc.c -std=c99 -I src -I deps -I deps/greatest -o $@
	@./$@

.PHONY: test
