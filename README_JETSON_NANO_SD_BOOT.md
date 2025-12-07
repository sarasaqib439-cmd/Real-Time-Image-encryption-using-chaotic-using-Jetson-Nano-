# Flashing Jetson Nano eMMC Module and Booting from SD Card

## Overview
Complete guide to flash a Jetson Nano eMMC module with NVIDIA L4T R32.7.6 and configure it to boot from SD card instead of eMMC.

## Prerequisites
- Jetson Nano Developer Kit with eMMC module
- Host PC running Ubuntu (tested on Ubuntu 20.04)
- Micro-USB cable for recovery mode
- HDMI display and serial console (optional but recommended for debugging)
- SD card (32GB+ recommended)

## Step 1: Download L4T Components

```bash
cd ~/nvidia
mkdir nvidia_sdk
cd nvidia_sdk

# Download L4T BSP (Board Support Package)
wget https://developer.nvidia.com/embedded/l4t/r32_release_v7.6/t210/jetson-210_linux_r32.7.6_aarch64.tbz2

# Download Sample Root Filesystem
wget https://developer.nvidia.com/embedded/l4t/r32_release_v7.6/t210/tegra_linux_sample-root-filesystem_r32.7.6_aarch64.tbz2

# Extract L4T BSP
tar xf Jetson-210_Linux_R32.7.6_aarch64.tbz2
```

## Step 2: Prepare Root Filesystem

```bash
cd Jetson-210_Linux_R32.7.6_aarch64/Linux_for_Tegra/rootfs

# Extract sample rootfs
sudo tar xpf ../../../Tegra_Linux_Sample-Root-Filesystem_R32.7.6_aarch64.tbz2

cd ..

# Apply NVIDIA binaries to rootfs
sudo ./apply_binaries.sh
```

## Step 3: Apply Critical Hardware Patch (PCN211181)

This patch fixes eMMC timing issues that cause boot hangs on certain Jetson Nano modules.

```bash
cd ~/nvidia/nvidia_sdk/Jetson-210_Linux_R32.7.6_aarch64/Linux_for_Tegra

# Download overlay patch
mkdir -p overlay
cd overlay
wget https://developer.nvidia.com/downloads/embedded/l4t/r32_release_v7.5/overlay_32.7.5_pcn211181.tbz2

# Extract and apply patch
tar xpf overlay_32.7.5_pcn211181.tbz2
cd ..
sudo cp -r overlay/* .

# Remove cached images
sudo rm -f bootloader/boot.img* bootloader/system.img*
```

## Step 4: Flash Jetson Nano eMMC

```bash
# Put Jetson Nano in Recovery Mode:
# 1. Power off the Jetson completely
# 2. Connect Micro-USB cable between Jetson and host PC
# 3. Short the FC REC (Force Recovery) pins with a jumper
# 4. Power on the Jetson
# 5. Remove the FC REC jumper

# Verify device is detected
lsusb | grep NVIDIA
# Should show: ID 0955:7f21 NVIDIA Corp. APX

# Flash the eMMC (use jetson-nano-devkit-emmc for developer kit with eMMC)
sudo ./flash.sh jetson-nano-devkit-emmc mmcblk0p1
```

**Note:** Flash will take 10-15 minutes. The Jetson will now boot from eMMC.

**Target Board Configurations:**
- `jetson-nano-emmc` - Standalone eMMC module (P3448 module only)
- `jetson-nano-devkit-emmc` - Developer Kit with eMMC
- `jetson-nano-devkit` - Developer Kit with SD card slot
- `jetson-nano-2gb-devkit` - 2GB Developer Kit

## Step 5: Enable SD Card Support

### 5.1: Compile SD Card Device Tree Overlay

On the host PC:

```bash
cd ~/nvidia/nvidia_sdk

# Clone device tree overlays repository
git clone https://github.com/Seeed-Studio/seeed-linux-dtoverlays.git
cd seeed-linux-dtoverlays

# Modify compatibility string for Jetson Nano
sed -i '17s#JETSON_COMPATIBLE#"nvidia,p3449-0000-b00+p3448-0002-b00", "nvidia,jetson-nano", "nvidia,tegra210"#' overlays/jetsonnano/jetson-sdmmc-overlay.dts

# Compile device tree overlay
make overlays/jetsonnano/jetson-sdmmc-overlay.dtbo
```

### 5.2: Install Overlay on Jetson

```bash
# Copy overlay to Jetson (replace IP with your Jetson's IP)
scp overlays/jetsonnano/jetson-sdmmc-overlay.dtbo nvidia@<JETSON_IP>:/tmp/

# SSH into Jetson
ssh nvidia@<JETSON_IP>
# Default credentials: username=nvidia, password=nvidia

# Install overlay
sudo cp /tmp/jetson-sdmmc-overlay.dtbo /boot/

# Configure hardware
sudo /opt/nvidia/jetson-io/config-by-hardware.py -n "reComputer sdmmc"

# Reboot to apply changes
sudo reboot
```

## Step 6: Prepare SD Card

After reboot, SSH back into the Jetson:

```bash
# Verify SD card is detected
lsblk
# Should show mmcblk1 (SD card) in addition to mmcblk0 (eMMC)

# Partition SD card
sudo fdisk /dev/mmcblk1
# In fdisk:
# - Press 'n' for new partition
# - Press 'p' for primary
# - Press '1' for partition number
# - Press Enter twice to accept defaults
# - Press 'w' to write changes

# Format SD card as ext4
sudo mkfs.ext4 -L ROOTFS_SD /dev/mmcblk1p1
```

## Step 7: Clone eMMC Rootfs to SD Card

```bash
# Create mount points
sudo mkdir -p /mnt/emmc /mnt/sd

# Mount both filesystems
sudo mount /dev/mmcblk0p1 /mnt/emmc
sudo mount /dev/mmcblk1p1 /mnt/sd

# Clone rootfs from eMMC to SD card (takes 5-10 minutes)
sudo rsync -aAXv --exclude={"/dev/*","/proc/*","/sys/*","/tmp/*","/run/*","/mnt/*","/media/*","/lost+found"} /mnt/emmc/ /mnt/sd/

# Unmount
sudo umount /mnt/emmc /mnt/sd
```

## Step 8: Configure Boot from SD Card

```bash
# Edit extlinux.conf to boot from SD card
sudo nano /boot/extlinux/extlinux.conf

# Change the APPEND line to use mmcblk1p1:
# Find the line that looks like:
#   APPEND ${cbootargs} root=/dev/mmcblk0p1 ...
# Change it to:
#   APPEND ${cbootargs} root=/dev/mmcblk1p1 rw rootwait rootfstype=ext4

# Save and reboot
sudo reboot
```

## Step 9: Verify Boot Source

After reboot:

```bash
# Check current root filesystem
mount | grep "on / "
# Should show: /dev/mmcblk1p1 on / type ext4 (rw,relatime,data=ordered)

# Check disk usage
df -h
# / should be mounted on /dev/mmcblk1p1 (SD card)
# /dev/mmcblk0p1 (eMMC) may be auto-mounted under /media

# Verify kernel command line
cat /proc/cmdline | grep -o "root=[^ ]*"
# Should show: root=/dev/mmcblk1p1
```

## Troubleshooting

### Boot Hangs at "bootconsole [uart8250] enabled"
- **Cause:** eMMC timing issues on affected hardware revisions
- **Solution:** Apply PCN211181 overlay patch (Step 3)
- **Symptoms:** System boots to bootloader, shows NVIDIA logo, but hangs before login prompt

### SD Card Not Detected (No mmcblk1)
- **Cause:** SD card controller disabled in device tree
- **Solution:** Install jetson-sdmmc-overlay.dtbo (Step 5)
- **Verify:** Check `dmesg | grep -i mmc` for SD card controller initialization
- **Common Issue:** Wrong Jetson configuration used during flash

### Wrong Root Device After Reboot
- **Verify:** Check `/boot/extlinux/extlinux.conf` has correct `root=/dev/mmcblk1p1`
- **Verify:** Run `cat /proc/cmdline` to see actual kernel boot parameters
- **Fix:** Edit `/boot/extlinux/extlinux.conf` and reboot

### Flash Command Fails with "Invalid target board"
- **Cause:** Using `.conf` extension or wrong configuration name
- **Solution:** Use configuration name without extension:
  - ✅ Correct: `sudo ./flash.sh jetson-nano-devkit-emmc mmcblk0p1`
  - ❌ Wrong: `sudo ./flash.sh jetson-nano-devkit-emmc.conf mmcblk0p1`

### Device Not in Recovery Mode
- **Verify:** `lsusb | grep NVIDIA` shows `ID 0955:7f21 NVIDIA Corp. APX`
- **Solution:** Re-enter recovery mode (power off, short FC REC, power on, remove short)

## Key Configuration Files

- **Boot Configuration:** `/boot/extlinux/extlinux.conf`
- **Device Tree Overlay:** `/boot/jetson-sdmmc-overlay.dtbo`
- **L4T Version:** Check with `cat /etc/nv_tegra_release`
- **Kernel Command Line:** `/proc/cmdline`

## Important Notes

### Boot Process
1. **TegraBoot** → Loads from eMMC boot partitions (mmcblk0boot0/1)
2. **cboot** → Bootloader reads extlinux.conf
3. **U-Boot** → Loads kernel and initrd
4. **Linux Kernel** → Mounts rootfs from device specified in `root=` parameter
5. **Systemd** → Initializes userspace

### Storage Layout
- **mmcblk0** = eMMC (14.7GB, non-removable)
  - mmcblk0p1 = Original rootfs partition
  - mmcblk0p2-p17 = Bootloader partitions
- **mmcblk1** = SD Card (size varies, removable)
  - mmcblk1p1 = Cloned rootfs partition

### Default Credentials
- **Username:** nvidia
- **Password:** nvidia

## Summary

You now have:
- ✅ Jetson Nano eMMC flashed with L4T R32.7.6
- ✅ PCN211181 patch applied (fixes boot hang issues)
- ✅ SD card controller enabled via device tree overlay
- ✅ Rootfs cloned from eMMC to SD card
- ✅ System booting from SD card (`/dev/mmcblk1p1`)
- ✅ eMMC still available as backup storage (`/dev/mmcblk0p1`)

The eMMC (14GB) remains accessible and can be mounted for additional storage or as a backup boot option.

## Additional Resources

- [NVIDIA Jetson Linux Developer Guide](https://docs.nvidia.com/jetson/archives/l4t-archived/l4t-3276/index.html)
- [L4T R32.7.6 Release Notes](https://developer.nvidia.com/embedded/linux-tegra-r3276)
- [PCN211181 Overlay Documentation](https://developer.nvidia.com/embedded/linux-tegra-r3275)
- [Seeed Studio Device Tree Overlays](https://github.com/Seeed-Studio/seeed-linux-dtoverlays)

## License

This guide is provided as-is for educational purposes. NVIDIA L4T and related components are subject to NVIDIA's licensing terms.
