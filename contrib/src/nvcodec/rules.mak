NVCODEC_HASH := 96a6db017b096ad48612890083464a7214902afa
NVCODEC_GITURL := https://git.videolan.org/git/ffmpeg/nv-codec-headers.git

ifndef HAVE_DARWIN_OS
PKGS += nvcodec
endif

$(TARBALLS)/nvcodec-$(NVCODEC_HASH).tar.xz:
	$(call download_git,$(NVCODEC_GITURL),,$(NVCODEC_HASH))

.sum-nvcodec: nvcodec-$(NVCODEC_HASH).tar.xz

nvcodec: nvcodec-$(NVCODEC_HASH).tar.xz .sum-nvcodec
	$(UNPACK)
	$(APPLY) $(SRC)/nvcodec/0001-allow-overriding-the-PREFIX-from-the-environment.patch
	$(MOVE)

.nvcodec: nvcodec-$(NVCODEC_HASH).tar.xz nvcodec
	cd nvcodec && PREFIX="$(PREFIX)" make install
	touch $@
