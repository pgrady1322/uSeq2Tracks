#!/usr/bin/env python3
"""
Quick test script to verify rapid mode configuration
"""

import yaml
import sys
import os

def test_rapid_mode_config():
    """Test that rapid mode configuration works correctly"""
    
    print("🧪 Testing uSeq2Tracks rapid mode configuration...")
    
    # Test 1: Load config file
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        print("✅ config.yaml loads successfully")
    except Exception as e:
        print(f"❌ Failed to load config.yaml: {e}")
        return False
    
    # Test 2: Check genome_id is present
    if 'genome_id' in config and config['genome_id']:
        print(f"✅ genome_id found: {config['genome_id']}")
    else:
        print("❌ genome_id missing or empty")
        return False
    
    # Test 3: Check rapid_mode parameter
    rapid_mode = config.get('rapid_mode', False)
    print(f"✅ rapid_mode setting: {rapid_mode}")
    
    # Test 4: Simulate output directory structure
    genome_id = config['genome_id'].strip().replace(" ", "_")
    base_outdir = config.get('outdir', 'results')
    genome_outdir = f"{base_outdir}/{genome_id}"
    
    print(f"📁 Output will be organized as: {genome_outdir}")
    
    if rapid_mode:
        rapid_outdir = f"{genome_outdir}/rapid"
        print(f"🚀 Rapid mode output directory: {rapid_outdir}")
        print("🚀 RAPID MODE: Will skip QC steps, generate essential tracks only")
    else:
        print("📊 STANDARD MODE: Will include full QC analysis")
    
    print("✅ All configuration tests passed!")
    return True

if __name__ == "__main__":
    if test_rapid_mode_config():
        print("\n🎉 uSeq2Tracks rapid mode configuration is ready!")
        sys.exit(0)
    else:
        print("\n❌ Configuration test failed")
        sys.exit(1)
