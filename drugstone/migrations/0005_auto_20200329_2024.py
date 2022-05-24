# Generated by Django 3.0.4 on 2020-03-29 20:24

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0004_protein_protein_description'),
    ]

    operations = [
        migrations.AlterField(
            model_name='edge',
            name='effect',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='drugstone.Effect'),
        ),
        migrations.AlterField(
            model_name='edge',
            name='protein_group',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='drugstone.ProteinGroup'),
        ),
    ]