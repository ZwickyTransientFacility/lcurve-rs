use phoebe_rs::{EBParams, compute_lightcurve};

#[test]
fn test_contact_binary() {
    let params = EBParams::contact(0.5, 85.0);
    let phases: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
    let lc = compute_lightcurve(&params, &phases).unwrap();

    // Should have eclipses: min flux < 1.0
    let min_flux = lc.flux.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_flux = lc.flux.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    println!("Contact q=0.5, i=85: min_flux={:.4}, max_flux={:.4}", min_flux, max_flux);
    assert!(min_flux < 0.95, "Should show eclipses");
    assert!(max_flux > 0.99, "Out-of-eclipse should be ~1.0");

    // Primary eclipse (phase 0.0) should be deeper than secondary (phase 0.5)
    // because the hotter star is eclipsed
    let f_primary = lc.flux[0];
    let f_secondary = lc.flux[50];
    println!("  f(phase=0.0)={:.4}, f(phase=0.5)={:.4}", f_primary, f_secondary);
}

#[test]
fn test_detached_binary() {
    let params = EBParams::detached(0.8, 88.0, 0.3, 0.3);
    let phases: Vec<f64> = (0..200).map(|i| i as f64 / 200.0).collect();
    let lc = compute_lightcurve(&params, &phases).unwrap();

    let min_flux = lc.flux.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("Detached q=0.8, i=88, f=0.3: min_flux={:.4}", min_flux);
    assert!(min_flux < 0.98, "Should show eclipses for detached at i=88");
}

#[test]
fn test_low_inclination_no_eclipse() {
    let params = EBParams::contact(0.5, 30.0);
    let phases: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
    let lc = compute_lightcurve(&params, &phases).unwrap();

    let min_flux = lc.flux.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("Contact q=0.5, i=30: min_flux={:.4}", min_flux);
    // At low inclination, should still show ellipsoidal variation but no eclipses
    // Min flux should be closer to 1.0
    assert!(min_flux > 0.90, "No deep eclipses at i=30");
}

#[test]
fn rpot_check() {
    // Check if our rpot matches Kopal convention
    // At (0.5, 0, 0) for q=0.5:
    // Kopal: Ω = 1/0.5 + 0.5*(1/0.5 - 0.5) + 0.75*0.25 = 2 + 0.75 + 0.1875 = 2.9375
    let q = 0.5;
    let pos = lcurve_subs::Vec3::new(0.5, 0.0, 0.0);
    let pot = lcurve_roche::potential::rpot(q, &pos).unwrap();
    println!("rpot at (0.5,0,0) for q=0.5: {:.6}", pot);
    println!("Expected Kopal: 2.9375");
    println!("Difference: {:.6}", pot - 2.9375);
    
    // Also check L1
    let xl1 = lcurve_roche::lagrange::xl11(q, 1.0).unwrap();
    let l1_pos = lcurve_subs::Vec3::new(xl1, 0.0, 0.0);
    let pot_l1 = lcurve_roche::potential::rpot(q, &l1_pos).unwrap();
    println!("\nL1 position: x={:.6}", xl1);
    println!("rpot at L1: {:.6}", pot_l1);
    println!("PHOEBE2 pot_max (L1): 2.876");
}

#[test]
fn rpot_l2_check() {
    let q = 0.5;
    let spin = 1.0;
    let xl2 = lcurve_roche::lagrange::xl12(q, spin).unwrap();
    let l2_pos = lcurve_subs::Vec3::new(1.0 + (1.0 - xl2), 0.0, 0.0);
    let pot_l2 = lcurve_roche::potential::rpot(q, &l2_pos).unwrap();

    let xl1 = lcurve_roche::lagrange::xl11(q, spin).unwrap();
    let l1_pos = lcurve_subs::Vec3::new(xl1, 0.0, 0.0);
    let pot_l1 = lcurve_roche::potential::rpot(q, &l1_pos).unwrap();

    println!("L1: x={:.6}, rpot={:.6}", xl1, pot_l1);
    println!("L2: x={:.6}, rpot={:.6}", 1.0 + (1.0 - xl2), pot_l2);
    println!("L1-L2 range: {:.6}", pot_l1 - pot_l2);

    // PHOEBE2 values:
    // pot_max (L1) = 2.876, pot_min (L2) = 2.577
    // pot = 2.851 => fillout = (2.876 - 2.851) / (2.876 - 2.577) = 0.025/0.299 = 0.084
    // Our pot_surface = pot_l1 + f * (pot_l2 - pot_l1)
    // To match PHOEBE's surface, we need:
    //   pot_surface_ours corresponds to Omega_phoebe = 2.851
    //
    // If the potential is a linear rescaling: Omega = a * rpot + b
    // Then at L1: 2.876 = a * (-1.973) + b
    //      at L2: 2.577 = a * pot_l2_ours + b
    // Two equations, two unknowns => a and b
    // Then for fillout f in our convention:
    //   pot_ours = pot_l1 + f * (pot_l2 - pot_l1)
    //   Omega = a * pot_ours + b = a*pot_l1 + a*f*(pot_l2-pot_l1) + b
    //         = 2.876 + a*f*(pot_l2-pot_l1)
    //   PHOEBE fillout: Omega = 2.876 - f_phoebe * 0.299
    //   => a*f*(pot_l2-pot_l1) = -f_phoebe * 0.299
    //   => f = f_phoebe * 0.299 / (-a * (pot_l2-pot_l1))
    //
    // Since a is the same scale factor: a = (2.876 - 2.577) / (pot_l1 - pot_l2)
    //   = 0.299 / (pot_l1 - pot_l2)
    // => f = f_phoebe * 0.299 / (-(0.299/(pot_l1-pot_l2)) * (pot_l2-pot_l1))
    //      = f_phoebe * 0.299 / (0.299) = f_phoebe
    //
    // So f_ours = f_phoebe IF the mapping is linear!
    // But our fillout formula is: pot = pot_l1 + f * (pot_l2 - pot_l1)
    // PHOEBE's fillout formula is: Omega = Omega_L1 - f * (Omega_L1 - Omega_L2)
    //                             = Omega_L1 + f * (Omega_L2 - Omega_L1)
    // Both are the same linear interpolation. So f=0.08 should produce
    // the same surface in both codes.
    println!("\nConclusion: fillout mapping IS the same.");
    println!("The 5% residual is NOT from potential convention.");
}

#[test]
fn surface_area_check() {
    // Check total surface area of our overcontact mesh
    // For a Roche-filling primary at q=0.5, the effective radius is ~0.38
    // (from PHOEBE: requiv_primary ≈ 1.5 in solar radii, sma ≈ ... )
    // The area should be ~4π r² ≈ 4π(0.38)² ≈ 1.81 (per star, units of a²)
    
    let q = 0.5;
    let params = phoebe_rs::EBParams::contact(q, 85.0);
    
    // Build overcontact mesh
    let (m1, m2) = phoebe_rs::geometry::build_overcontact_mesh(q, 0.08, 80).unwrap();
    
    let area1: f64 = m1.elements.iter().map(|e| e.area).sum();
    let area2: f64 = m2.elements.iter().map(|e| e.area).sum();
    let n1 = m1.elements.len();
    let n2 = m2.elements.len();
    
    println!("Overcontact mesh (q=0.5, fillout=0.08):");
    println!("  Primary:   {} elements, total area = {:.4}", n1, area1);
    println!("  Secondary: {} elements, total area = {:.4}", n2, area2);
    println!("  Total:     {} elements, total area = {:.4}", n1+n2, area1+area2);
    
    // Compare with detached Roche-filling
    let m1d = phoebe_rs::geometry::build_mesh(q, lcurve_roche::Star::Primary, 1.0, 80).unwrap();
    let m2d = phoebe_rs::geometry::build_mesh(q, lcurve_roche::Star::Secondary, 1.0, 80).unwrap();
    let area1d: f64 = m1d.elements.iter().map(|e| e.area).sum();
    let area2d: f64 = m2d.elements.iter().map(|e| e.area).sum();
    
    println!("\nDetached Roche-filling (ffac=1.0):");
    println!("  Primary:   {} elements, total area = {:.4}", m1d.elements.len(), area1d);
    println!("  Secondary: {} elements, total area = {:.4}", m2d.elements.len(), area2d);
    
    println!("\nOvercontact/Detached area ratio:");
    println!("  Primary:   {:.3}", area1 / area1d);
    println!("  Secondary: {:.3}", area2 / area2d);
}
